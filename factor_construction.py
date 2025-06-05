import numpy as np
import scipy as sp
import random
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from time import time
import itertools
import multiprocessing as mp
import sys
from contextlib import redirect_stdout

from configs import *

# Avoid numpy multithreading conflicts
import os
os.environ['OMP_NUM_THREADS'] = '1'
os.environ['OPENBLAS_NUM_THREADS'] = '1' 
os.environ['MKL_NUM_THREADS'] = '1'
os.environ['VECLIB_MAXIMUM_THREADS'] = '1'
os.environ['NUMEXPR_NUM_THREADS'] = '1'


def match_surveys_with_subgroups(surveys_2, search_groups= [[7,64,73,55], [43,46,51,50]]):
    """Searches through surveys to find all those that match the given subgroups."""
    match_row_idx = []
    for i, row in enumerate(surveys_2):
        subgroup_1_match = all(item in row[0:5] for item in search_groups[0]) or all(item in row[5:10] for item in search_groups[0])
        subgroup_2_match = all(item in row[0:5] for item in search_groups[1]) or all(item in row[5:10] for item in search_groups[1])

        if subgroup_1_match and subgroup_2_match:
            match_row_idx.append(i)

    # Get the average scores for these rows
    match_surveys = surveys_2[match_row_idx, :].astype(int)

    # Remove duplicate rows with different orderings
    unique_surveys = []
    for survey in match_surveys:
        # Order the groups
        ordered_survey = np.sort(survey[0:5]).tolist() + np.sort(survey[5:10]).tolist()
        if ordered_survey not in unique_surveys:
            unique_surveys.append(ordered_survey)

    return match_row_idx, np.array(unique_surveys)

def get_groups_from_keep(keep, n_per_group):
    return [keep[i:i + n_per_group] for i in range(0, len(keep), n_per_group)]

def get_best_survey_info(surveys, avg_score, n_per_group):
    """Sorts surveys by average score and returns the best survey, its groups, and the sorting indices."""
    # Sort on average score
    sorted_inds = np.flip(np.argsort(avg_score))

    # Get the items to keep and the groups
    keep = (surveys[sorted_inds[0],:]).astype(int).tolist()
    groups = get_groups_from_keep(keep, n_per_group)

    return keep, groups, sorted_inds

def iterative_factor_construction(ansbyq, seed_indices, n_per_group=None, n_subj_drop=0, drop_samples=0):
    """
    Iteratively construct factors by expanding seed groups based on correlation dominance.

    Arguments:
        seed_indices: List of sublists, each sub-list contains question indices to seed a cluster
        n_per_group: Stop when each group reaches this size. If None, continues until all items assigned.
    
    Returns:
        groups: List of sublists, each sublist is a set of indices for a group (factor)
        fit_history: Dictionary of fit scores for items added at each step
    """
    n_questions = ansbyq.shape[0]
    n_subjects = ansbyq.shape[1]
  
    # Calculate correlation matrix, averaging over samples if requested
    corr_matrix = np.zeros((n_questions, n_questions))
    if drop_samples > 0:
        for i in range(drop_samples):
            # Randomly drop subjects
            keep = np.random.choice(range(n_subjects), size=n_subjects - n_subj_drop, replace=False)
            sample = ansbyq[:, keep]
            sample_corr = np.corrcoef(sample)
            corr_matrix += sample_corr
    
        # Normalize by number of samples
        corr_matrix /= drop_samples
    else:
        # Use full correlation matrix
        corr_matrix = np.corrcoef(ansbyq)

    # Initialize groups with seeds
    groups = [list(seed) for seed in seed_indices]
    assigned_items = set(item for group in groups for item in group)
    available_items = set(range(n_questions)) - assigned_items
    
    # Set stopping criterion
    if n_per_group is None:
        n_per_group = n_questions // len(groups)  # Distribute evenly
    
    # Initialize fit history
    fit_history = []
    
    # Continue until stopping criterion met
    while any(len(group) < n_per_group for group in groups) and available_items:
        
        # Cycle through groups
        for group_idx, group in enumerate(groups):
            if len(group) >= n_per_group or not available_items:
                continue
                
            best_item = None
            best_fit_score = -1
            
            # Evaluate each available item for this group
            for item in available_items:
                fit_scores_all_groups = []
                
                # Calculate average correlation with each group
                for other_group in groups:
                    if len(other_group) > 0:
                        group_corrs = [abs(corr_matrix[item, member]) for member in other_group]
                        avg_corr = np.mean(group_corrs)
                        fit_scores_all_groups.append(avg_corr)
                    else:
                        fit_scores_all_groups.append(0)
                
                # Calculate dominance ratio 
                target_group_corr = fit_scores_all_groups[group_idx]
                total_corr = np.sum(fit_scores_all_groups)
                
                # Avoid division by zero
                if total_corr > 0:
                    dominance_ratio = target_group_corr / total_corr
                else:
                    dominance_ratio = 0
                
                # Track best item for this group
                if dominance_ratio > best_fit_score:
                    best_fit_score = dominance_ratio
                    best_item = item
            
            # Add best item to group
            if best_item is not None:
                groups[group_idx].append(best_item)
                available_items.remove(best_item)
                assigned_items.add(best_item)
                
                fit_history.append({
                    'group': group_idx,
                    'item': best_item,
                    'fit_score': best_fit_score
                })
    
    return groups, fit_history


def assess_group_quality(ansbyq, groups, fit_history, skip_pca=True):
    """Assess quality of groups using aggregate fit scores and PCA variance explained."""
    
    n_groups = len(groups)
    
    # Calculate aggregate fit scores for each group
    group_fit_scores = [[] for _ in range(n_groups)]
    
    # Collect fit scores by group
    for entry in fit_history:
        group_idx = entry['group']
        fit_score = entry['fit_score']
        group_fit_scores[group_idx].append(fit_score)
    
    # Calculate mean fit score for each group
    aggregate_fit_scores = []
    for group_scores in group_fit_scores:
        if len(group_scores) > 0:
            aggregate_fit_scores.append(np.mean(group_scores))
        else:
            aggregate_fit_scores.append(0.0)  # Handle empty groups
    
    # Subset ansbyq to questions that ended up in groups
    all_group_items = [item for group in groups for item in group]
    
    if len(all_group_items) == 0:
        return aggregate_fit_scores, 0.0
    
    subset_data = ansbyq[all_group_items, :]  # subset rows (questions)
    
    # Perform PCA (transpose for sklearn: samples as rows, features as columns)
    pca_data = subset_data.T
    
    # Determine number of components to extract
    max_components = min(n_groups, pca_data.shape[1], pca_data.shape[0])
    
    if max_components <= 0:
        return aggregate_fit_scores, 0.0
    
    # For exhaustive checking, too expensive
    if skip_pca:
        cumulative_variance_explained = np.nan
    else:
        pca = PCA(n_components=max_components)
        pca.fit(pca_data)
        
        # Calculate cumulative variance explained by first n_groups components
        n_components_to_sum = min(n_groups, len(pca.explained_variance_ratio_))
        cumulative_variance_explained = np.sum(pca.explained_variance_ratio_[:n_components_to_sum])
    
    return aggregate_fit_scores, cumulative_variance_explained

def assess_factor_stability(ansbyq, keep, subjects_to_drop_list, n_samples=100, n_components=None):
    """Assess stability of factor structure through subsampling."""
    
    # Subset to kept questions
    subset_data = ansbyq[keep, :]
    n_questions, n_subjects = subset_data.shape
    
    # Fit reference PCA on full subset
    reference_data = subset_data.T  # transpose for sklearn
    
    if n_components is None:
        n_components = min(n_questions, n_subjects, 10)  # cap at 10 for practical purposes
    
    # Handle edge case where we can't fit the requested components
    max_possible_components = min(reference_data.shape[0], reference_data.shape[1])
    n_components = min(n_components, max_possible_components)
    
    if n_components <= 0:
        return {'error': 'Insufficient data for PCA'}
    
    reference_pca = PCA(n_components=n_components)
    reference_pca.fit(reference_data)
    reference_loadings = reference_pca.components_  # shape: (n_components, n_questions)
    
    stability_results = {
        'subjects_dropped': [],
        'mean_congruence': [],
        'std_congruence': [],
        'sample_congruences': [],
        'reference_variance_explained': reference_pca.explained_variance_ratio_
    }
    
    for n_drop in subjects_to_drop_list:
        if n_drop >= n_subjects:
            # Skip if we'd drop all subjects
            continue
            
        remaining_subjects = n_subjects - n_drop
        
        if remaining_subjects < n_components:
            # Skip if insufficient subjects for PCA
            continue
        
        sample_congruences = []
        
        for sample_idx in range(n_samples):
            # Random sample of subjects to keep
            keep_subjects = random.sample(range(n_subjects), remaining_subjects)
            
            # Subset data
            sample_data = subset_data[:, keep_subjects].T  # transpose for sklearn
            
            # Fit PCA on sample
            try:
                sample_pca = PCA(n_components=n_components)
                sample_pca.fit(sample_data)
                sample_loadings = sample_pca.components_
                
                # Calculate Tucker's congruence coefficient for each component
                component_congruences = []
                for comp_idx in range(n_components):
                    ref_comp = reference_loadings[comp_idx, :]
                    sample_comp = sample_loadings[comp_idx, :]
                    
                    # Handle sign indeterminacy by taking absolute value of correlation
                    congruence = abs(np.corrcoef(ref_comp, sample_comp)[0, 1])
                    
                    # Handle NaN case (constant components)
                    if np.isnan(congruence):
                        congruence = 0.0
                        
                    component_congruences.append(congruence)
                
                # Average congruence across components
                mean_sample_congruence = np.mean(component_congruences)
                sample_congruences.append(mean_sample_congruence)
                
            except Exception as e:
                # PCA failed for this sample
                sample_congruences.append(0.0)
        
        # Store results for this drop level
        stability_results['subjects_dropped'].append(n_drop)
        stability_results['mean_congruence'].append(np.mean(sample_congruences))
        stability_results['std_congruence'].append(np.std(sample_congruences))
        stability_results['sample_congruences'].append(sample_congruences)
    
    # Convert lists to numpy arrays for easier handling
    stability_results['subjects_dropped'] = np.array(stability_results['subjects_dropped'])
    stability_results['mean_congruence'] = np.array(stability_results['mean_congruence'])
    stability_results['std_congruence'] = np.array(stability_results['std_congruence'])
    stability_results['sample_congruences'] = np.array(stability_results['sample_congruences'])

    return stability_results

def evaluate_two_factor_surveys(ansbyq, combs=None, n_per_group=5, verbose=0, n_subj_drop=0, drop_samples=0):
    """
    Evaluate all 2-group seeds, returning all generated surveys, along with info from the best.
    Note that it would take a long time to evaluate all 3-group seeds.
    n_possible_3_group_seeds = sp.special.comb(100,3) = 161700
    """
    n_groups = 2

    # Default to exhaustive enumeration
    if combs is None:
        combs = itertools.combinations(np.arange(0,100,1), n_groups)
        n_samples = int(sp.special.comb(100, n_groups))
    else:
        n_samples = len(combs)

    # Initialize data arrays
    surveys     = np.zeros((n_samples, n_groups* n_per_group), dtype=int)
    scores      = np.zeros((n_samples, n_groups))
    avg_score   = np.zeros((n_samples))
    varexp      = np.zeros((n_samples))

    # Evaluate all 2-group seeds
    start_time = time()
    for i, comb in enumerate(combs):
        
        # Set the seeds for the iterative construction
        seeds = [[i] for i in comb]

        # Perform iterative factor construction, assess survey quality
        groups, fit_history = iterative_factor_construction(ansbyq, seeds, n_per_group=n_per_group, n_subj_drop=n_subj_drop, drop_samples=drop_samples)
        aggregate_fit_scores, cumulative_variance_explained = assess_group_quality(ansbyq, groups, fit_history)

        # Store results
        surveys[i, :]   = np.array([i for grp in groups for i in grp], dtype=int)
        scores[i , :]   = aggregate_fit_scores
        avg_score[i]    = np.mean(aggregate_fit_scores)
        varexp[i]       = cumulative_variance_explained

        # Print progress if verbose
        if verbose > 0:
            curtime = time() - start_time
            estimated_end_time = (curtime * (n_samples - i) / (i + 1))/ 60  # in minutes
            print(f"Iteration {i+1}/{n_samples} completed after {curtime:.2f} seconds. Estimated time remaining: {estimated_end_time:.2f} min.")

    return surveys, avg_score, scores, varexp

def _evaluate_single_sample(args):
    """Worker function for parallel evaluation of all 2-factor surveys generated using an excluded-column subset of ansbyq."""
    ansbyq, drop_size, combs, sample_idx, random_seed = args
    
    # Set random seed for reproducibility
    np.random.seed(random_seed)
    
    # Randomly drop subjects
    n_subjects = ansbyq.shape[1]
    keep = np.random.choice(range(n_subjects), size=n_subjects - drop_size, replace=False)
    sub_survey = ansbyq[:, keep]
    
    # Evaluate the survey with the dropped subjects
    surveys, avg_score, _, _ = evaluate_two_factor_surveys(sub_survey, combs=combs, verbose=0)
    
    return surveys, avg_score, sample_idx

def evaluate_surveys_over_dropped_subjects(ansbyq, combs=None, drop_subj_list=[10], n_samples=20, n_processes=None, random_seed=42):
    """Drops random subject ids and evaluates 2- survey scores, aggregates results."""
    # Parameters: Don't change n_groups, set n_per_group in 2-6
    n_groups = 2
    n_per_group = 5

    # Default to exhaustive enumeration
    if combs is None:
        combs = [comb for comb in itertools.combinations(np.arange(0, 100, 1), n_groups)]

    # Set number of processes
    if n_processes is None:
        n_processes = mp.cpu_count()

    # Initialize storage for all aggregate scores
    avg_scores_aggregated = np.zeros((len(drop_subj_list), n_samples, len(combs)))

    # Initialize storage for all resulting surveys
    surveys_aggregated = np.zeros((len(drop_subj_list), n_samples, len(combs), n_groups*n_per_group), dtype=int)

    # Loop over drop counts
    for i, drop_size in enumerate(drop_subj_list):
        print(f"Processing drop size {drop_size} with {n_samples} samples using {n_processes} processes...")
        
        # Prepare arguments for all samples for this drop_size
        args_list = []
        for sample_idx in range(n_samples):
            # Generate unique random seed for each sample
            sample_seed = random_seed + i * n_samples + sample_idx
            args_list.append((ansbyq, drop_size, combs, sample_idx, sample_seed))
        
        # Execute samples in parallel
        with mp.Pool(n_processes) as pool:
            results = pool.map(_evaluate_single_sample, args_list)
        
        # Collect results
        for surveys, avg_score, sample_idx in results:
            # Surveys are ordered by seeds, so save groups resulting from each seed
            surveys_aggregated[i, sample_idx, :, :] = surveys
            # Store the average scores for each survey
            avg_scores_aggregated[i, sample_idx, :] = avg_score
        
        print(f"Completed drop size {drop_size}")

    return avg_scores_aggregated, surveys_aggregated



def evaluate_three_factor_surveys(ansbyq, keep, groups_in):
    """Evaluate all third factors from greedy 2-group solution"""
    remaining = np.setdiff1d(range(100), keep)
    n_samples = len(remaining)

    n_groups = 3
    n_per_group = 5
    surveys = np.zeros((n_samples, n_groups* n_per_group))
    scores  = np.zeros((n_samples, n_groups))
    avg_score = np.zeros((n_samples))
    varexp  = np.zeros((n_samples))

    start_time = time()
    for i, seed in enumerate(remaining):
        seeds = groups_in + [[seed]]

        groups, fit_history = iterative_factor_construction(ansbyq, seeds, n_per_group=n_per_group)
        aggregate_fit_scores, cumulative_variance_explained = assess_group_quality(ansbyq, groups, fit_history)

        surveys[i, :] = np.array([i for grp in groups for i in grp])
        avg_score[i] = np.mean(aggregate_fit_scores)
        scores[i , :] = aggregate_fit_scores
        varexp[i] = cumulative_variance_explained

        curtime = time() - start_time
        estimated_end_time = (curtime * (n_samples - i) / (i + 1))/ 60  # in minutes

        print(f"Iteration {i+1}/{n_samples} completed after {curtime:.2f} seconds. Estimated time remaining: {estimated_end_time:.2f} min.")

    sorted_inds = np.flip(np.argsort(avg_score))

    # Convert back to sublists of indices
    keep = (surveys[sorted_inds[0],:]).astype(int).tolist()
    groups = [keep[i:i + n_per_group] for i in range(0, len(keep), n_per_group)]

    return keep, groups, surveys, avg_score, scores, varexp


def analyze_seed_stability(surveys_aggregated, combs, drop_size_idx=0, n_per_group=5):
    """Analyze the stability of surveys across repetitions for each seed."""

    # Get the parameters of the aggregated survey data
    n_drop_sizes, n_samples, n_seeds, total_items = surveys_aggregated.shape
    n_groups = total_items // n_per_group
    
    stability_details = {}
    jaccard_avg = np.zeros(n_seeds)
    
    for seed_idx, seed in enumerate(combs):
        # Get all survey repetitions for this seed at specified drop size
        # Shape: (n_samples, total_items)
        seed_surveys = surveys_aggregated[drop_size_idx, :, seed_idx, :]
        
        # Reshape to separate groups: (n_samples, n_groups, n_per_group)
        seed_surveys_grouped = seed_surveys.reshape(n_samples, n_groups, n_per_group)
        
        # Count frequency of each item across all repetitions
        all_items = seed_surveys.flatten()
        unique_items, counts = np.unique(all_items, return_counts=True)
        
        # Calculate Jaccard index for every pair of repetitions
        jaccard_scores = []
        for i in range(n_samples):
            for j in range(i + 1, n_samples):
                set_i = set(seed_surveys[i])
                set_j = set(seed_surveys[j])
                intersection = len(set_i & set_j)
                union = len(set_i | set_j)
                jaccard = intersection / union if union > 0 else 0
                jaccard_scores.append(jaccard)
        
        # For each group, determine item frequencies across repetitions
        group_item_frequencies = []
        for group_idx in range(n_groups):

            # For each group, we'll check each item and the number of repetitions it appeared in
            group_item_counts = {}
            for rep_idx in range(n_samples):
                for item in seed_surveys_grouped[rep_idx, group_idx, :]:
                    group_item_counts[item] = group_item_counts.get(item, 0) + 1
            
            # Calculate what percentage of repetitions each item appeared in this group
            item_frequencies = {item: count/n_samples*100 for item, count in group_item_counts.items()}
            group_item_frequencies.append(sorted(item_frequencies.items(), key=lambda x: x[1], reverse=True))
        
        # Store results for this seed
        stability_details[seed_idx] = {
            'seed': seed,
            'jaccard_avg': np.mean(jaccard_scores) if jaccard_scores else 0,
            'jaccard_std': np.std(jaccard_scores) if jaccard_scores else 0,
            'group_item_frequencies': group_item_frequencies,
            'n_unique_items_seen': len(unique_items)
        }

        # Store average Jaccard index for this seed
        jaccard_avg[seed_idx] = stability_details[seed_idx]['jaccard_avg']

    return jaccard_avg, stability_details


def print_seed_info(seed_idx, qs, score_info, stability_details, top_n=5):
    """Print detailed information about a specific seed's stability and item frequencies."""

    print(f"Seed details for seed {stability_details[seed_idx]['seed']}:")
    print(f"Weighted Score  : {score_info['weighted_score'][seed_idx]:.1f}")
    print(f"Stability Score : {score_info['prctile_stability'][seed_idx]:.1f}")
    print(f"Fit Score       : {score_info['prctile_scores'][seed_idx]:.1f}")
    print(f"Top {top_n} most frequent items in each group:")
    for group_idx, group_freq in enumerate(stability_details[seed_idx]['group_item_frequencies']):
        print(f"  Group {group_idx + 1}:")
        for item, freq in group_freq[:top_n]:
            item = int(item)  # Ensure item is an integer index
            if item < len(qs):
                question_text = qs[item]
                print(f"    Item {item}: {freq:.1f}% - \"{question_text}\"")
            else:
                print(f"    Item {item}: {freq:.1f}% - [Question text not available]")
    print(' ')



def plot_stability_analysis(ansbyq, keep, n_components):
    """Performs subsetting analysis by dropping subjects and assessing factor stability."""
    drop_sizes = np.arange(0,20,1)
    stability_results = assess_factor_stability(ansbyq, keep, drop_sizes, n_components=n_components)

    # Plot the stability results
    plt.figure()
    plt.errorbar(drop_sizes, stability_results['mean_congruence'], yerr=2*stability_results['std_congruence'], fmt='o')
    plt.xticks(drop_sizes)
    plt.xlabel('Number of Subjects Dropped')
    plt.ylabel('Mean Congruence Coefficient')
    plt.title(f'Stability of {n_components}-Factor Structure Given Dropped Subjects')
    plt.grid()
    plt.legend(['Mean Congruence ±2 SD'])
    plt.tight_layout()

def print_group_questions(qs, groups):
    """Prints questions in each group."""
    for i, group in enumerate(groups):
        print(' ')
        print("Group questions: ", group)
        for q in [qs[i] for i in group]: 
            print(q)

def plot_subsurvey(ansbyq, qs, keep, n_groups):
    """Plot PCA of the selected questions in the survey."""

    newbyq = ansbyq[keep, :]
    newqs = [qs[i][0:20] for i in keep]

    pca = PCA()
    pca.fit(newbyq.T)

    plt.figure(figsize=(6, 5))
    for i in range(n_groups):
        plt.plot(newqs, pca.components_[i, :], '-o', label=f'Factor {i+1}')
    plt.axhline(0, color='gray', linestyle='--')
    plt.xticks(rotation=90)
    plt.xlabel('Questions')
    plt.ylabel('PCA Component Loadings')
    plt.title('PCA Loadings for Selected Questions')
    plt.legend()
    plt.tight_layout()

    plt.figure()
    plt.plot(pca.explained_variance_ratio_, 'o-')
    plt.xlabel('Principal Component Index')
    plt.ylabel('Explained Variance Ratio')
    plt.title('Explained Variance Ratio by Principal Component')
    plt.tight_layout()

def plot_subsurvey_construct_domain_heatmaps(qstats, qcoords, groups):
    """Wrapper for plotting construct-domain item grids given groups."""
    for group in groups:
        qstats_f = qstats.copy()
        to_nan = np.setdiff1d(range(0,100), group)
        qstats_f.iloc[to_nan] = np.nan

        mean_heatmap = create_construct_domain_item_heatmap(qstats_f['mean'].values, qcoords)
        plot_construct_domain_item_heatmap(mean_heatmap, vmin=1, vmax=5, title='Mean Values')

def plot_stability_vs_average_score(jaccard_avg, scores_avg, combs):
    """Plot stability (Jaccard index) vs average score as a scatter plot."""
    
    # Create scatter plot
    plt.figure(figsize=(10, 8))
    scatter = plt.scatter(scores_avg, jaccard_avg, alpha=0.6)
    
    # Add labels and title
    plt.xlabel('Mean Average Score')
    plt.ylabel('Stability (Mean Jaccard Index)')
    plt.title('Survey Stability vs Average Score')
    plt.grid(True, alpha=0.3)
    
    # Add correlation coefficient
    correlation = np.corrcoef(scores_avg, jaccard_avg)[0, 1]
    plt.text(0.05, 0.95, f'Correlation: {correlation:.3f}', 
             transform=plt.gca().transAxes, fontsize=12)
    
    # Optionally label some extreme points
    # Find highest scoring and most stable seeds
    top_score_idx   = np.argmax(scores_avg)
    most_stable_idx = np.argmax(jaccard_avg)
    
    plt.annotate(f'Best Score: {combs[top_score_idx]}', 
                xy=(scores_avg[top_score_idx], jaccard_avg[top_score_idx]),
                xytext=(10, 10), textcoords='offset points')
    
    plt.annotate(f'Most Stable: {combs[most_stable_idx]}', 
                xy=(scores_avg[most_stable_idx], jaccard_avg[most_stable_idx]),
                xytext=(10, -20), textcoords='offset points')
    
    plt.tight_layout()
    plt.show()
    
    # Print some summary statistics
    print(f"Correlation between stability and average score: {correlation:.3f}")
    print(f"Best scoring seed: {combs[top_score_idx]} (Score: {scores_avg[top_score_idx]:.3f}, Stability: {jaccard_avg[top_score_idx]:.3f})")
    print(f"Most stable seed : {combs[most_stable_idx]} (Score: {scores_avg[most_stable_idx]:.3f}, Stability: {jaccard_avg[most_stable_idx]:.3f})")


def get_weighted_survey_quality(jaccard_avg, scores_avg, weights=[0.5, 0.5]):
    """Calculate a weighted survey quality score based on stability and average score."""    
    # Convert jaccard_avg and scores_avg to percentiles
    prctile_stability = np.argsort(np.argsort(jaccard_avg)) / len(jaccard_avg) * 100
    prctile_scores    = np.argsort(np.argsort(scores_avg )) / len(scores_avg ) * 100

    # Calculate weighted score
    weighted_score = weights[0] * prctile_stability + weights[1] * prctile_scores
    return weighted_score, prctile_stability, prctile_scores


def get_consensus_items(best_indices, stability_details, top_n=5):
    """Creates a list of surveys based on the best_indices provided into the stability_detials dictionary."""
    n_groups = 2
    consensus_items = np.zeros((len(best_indices), top_n*n_groups), dtype=int)
    for j, indices in enumerate(best_indices):
        for group_idx, group_freq in enumerate(stability_details[indices]['group_item_frequencies']):
            for i, (item, freq) in enumerate(group_freq[:top_n]):
                consensus_items[j, i+ group_idx*top_n] = item
    
    return consensus_items


def run_factor_construction():
    # Drop subjects before or after construction of sub-surveys?
    drop_post_construction = False

    # Read subject data
    sids, data = read_response_data(subj_data_dir)

    # Chunk question data, get basic stats
    qs, ansbyq, qstats, qcdfs, qpdfs, inds = get_qdata(sids, data, sort=False, trim=True, flip=True)

    # List the seeds for the 2-factor surveys
    combs = [comb for comb in itertools.combinations(np.arange(0, 100, 1), 2)]

    if drop_post_construction:
        # Evaluate all 2-factor surveys in terms of stability and fit scores
        avg_scores_aggregated, surveys_aggregated = evaluate_surveys_over_dropped_subjects(ansbyq, n_processes=mp.cpu_count(), n_samples=25)
    else:
        # Evaluate all 2-factor surveys in terms of stability and fit scores
        surveys, avg_score, _, _ = evaluate_two_factor_surveys(ansbyq, combs=combs, n_per_group=5, verbose=1, n_subj_drop=0, drop_samples=0)

        # Reshape avg_score to match avg_scores_aggregated shape
        n_samples = 1  # Since we didn't drop subjects, we only have one sample
        avg_scores_aggregated = np.zeros((1, n_samples, len(combs)))
        avg_scores_aggregated[0, 0, :] = avg_score

        # Reshape surveys to match surveys_aggregated shape
        surveys_aggregated = np.zeros((1, n_samples, len(combs), 10))  # 10 items per group
        surveys_aggregated[0, 0, :, :] = surveys

        # Set weights to use only average score 
        weights = [0, 1]

    # Analyze stability for the first drop size
    jaccard_avg, stability_details = analyze_seed_stability(surveys_aggregated, combs, drop_size_idx=0, n_per_group=5)

    # Get average scores over repetitions using different subjects
    scores_avg = np.mean(avg_scores_aggregated[0,:,:], axis=0)

    if drop_post_construction:
        # Plot stability vs average score
        plot_stability_vs_average_score(jaccard_avg, scores_avg, combs)

    # Get score and stability percentiles, and weighted sum of these
    weighted_score, prctile_stability, prctile_scores = get_weighted_survey_quality(jaccard_avg, scores_avg, weights=weights)

    # Score information
    score_info = {
        'weighted_score': weighted_score,
        'prctile_stability': prctile_stability,
        'prctile_scores': prctile_scores,
        }

    # Get best seeds based on weighted score
    n_top_surveys = 5
    sorted_inds = np.flip(np.argsort(weighted_score))
    best_indices  = sorted_inds[:n_top_surveys]

    # Print info from the best seeds
    for seed_idx in best_indices:
        print_seed_info(seed_idx, qs, score_info, stability_details, top_n=5)

    # Get consensus groups from the best seeds
    consensus_items = get_consensus_items(best_indices, stability_details, top_n=5)

    # Perform PCA for each of the consensus surveys
    for i, items in enumerate(consensus_items):
        print(f"Consensus Survey {i+1}: {items}")
        plot_subsurvey(ansbyq, qs, items, n_groups=2)


    # keep_3, groups_3, surveys_3, avg_score_3, scores_3, varexp_3 = evaluate_three_factor_surveys(ansbyq, keep_2, groups_2)

    # plot_subsurvey(ansbyq, qs, keep_3, groups_3)
    # print_group_questions(qs, groups_3)

    # plot_stability_analysis(ansbyq, keep_2, n_components=2)
    # plot_stability_analysis(ansbyq, keep_3, n_components=3)

    # n_surveys = int(sp.special.comb(100, 2))

    # sorted_inds = np.argsort(avg_score_2)

    # # assessment set
    # assessment_inds = sorted_inds[np.arange(0, n_surveys, 5)]

    # # Asessment seeds
    # combs = []
    # for i, comb in enumerate(itertools.combinations(np.arange(0,100,1), 2)):
    #     if i in assessment_inds:
    #         combs.append(comb)


    # avg_avg_scores = np.mean(avg_scores_aggregated, axis=1)
    # std_avg_scores = np.std(avg_scores_aggregated, axis=1)

    # assessment_avg_scores_sort_idx = np.argsort(avg_score_2[np.sort(assessment_inds)])

    # plt.figure()
    # plt.plot(avg_avg_scores[0,assessment_avg_scores_sort_idx], label='5 Subjects Dropped')

    # corrs=[]
    # for i, drop_size in enumerate(drop_subj_list):
    #     corrs.append(sp.stats.spearmanr(avg_score_2[assessment_inds],avg_avg_scores[i, assessment_avg_scores_sort_idx])[0])