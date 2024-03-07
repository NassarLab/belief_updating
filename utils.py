import numpy as np

def get_cdf(data, type = 'cdf'):
    vals = np.unique(data)
    cnts = np.array([sum(item == data) for item in vals])
    csum = np.cumsum(cnts)
    cdf  = csum/csum[-1]
    pdf  = cnts/sum(cnts)

    df = cdf if type == 'cdf' else pdf

    return vals, df
    

plt.figure()
vals, cdf = get_cdf(dem_1['Age'])
plt.plot(vals, cdf)
vals, cdf = get_cdf(dem_2['Age'])
plt.plot(vals, cdf)

plt.figure()
vals, cdf = get_cdf(dem_1['Time taken'])
plt.plot(vals, cdf)
vals, cdf = get_cdf(dem_2['Time taken'])
plt.plot(vals, cdf)
vals, cdf = get_cdf(dem_1['Time taken']*(94/150))
plt.plot(vals, cdf)




plt.figure()
vals, cdf = get_cdf(dem_1['Sex'])
plt.plot(vals, cdf)
vals, cdf = get_cdf(dem_2['Sex'])
plt.plot(vals, cdf)


plt.figure()
vals, cdf = get_cdf(dem_1['Total approvals'])
plt.plot(vals, cdf)
vals, cdf = get_cdf(dem_2['Total approvals'])
plt.plot(vals, cdf)


plt.figure()
vals, cdf = get_cdf(dem_1['Ethnicity simplified'])
plt.plot(vals, cdf)
vals, cdf = get_cdf(dem_2['Ethnicity simplified'])
plt.plot(vals, cdf)


plt.figure();
plt.plot(dem_1.loc[dem_1['Sex'] == 'Female']['Age'],'o')
plt.plot(dem_1.loc[dem_1['Sex'] == 'Male']['Age'],'o')

plt.figure();
plt.plot(dem_2.loc[dem_2['Sex'] == 'Female']['Age'],'o')
plt.plot(dem_2.loc[dem_2['Sex'] == 'Male']['Age'],'o')



plt.figure(); 
plt.plot(np.arange(94), qstats_1.loc[match_pairs[:,1]+1,:]['mean'],'o')
plt.plot(np.arange(94), qstats_2['mean'],'o')



plt.figure(); 
plt.plot(np.arange(94), qstats_1.loc[match_pairs[:,1]+1,:]['std'],'o')
plt.plot(np.arange(94), qstats_2['std'],'o')