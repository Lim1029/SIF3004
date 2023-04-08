import matplotlib.pyplot as plt
import pandas as pd

# de-boosting equation:

def deboost(s_obs,snr):
    s_true = s_obs / (1+0.2*(snr/5)**-2.3)
    return s_true

catalog_path = '/home/mingkang/Desktop/bootes/sources.dat'
data = pd.read_table(catalog_path,sep='\s+',index_col='obj_id')
data['deboost'] = deboost(data['flux'],data['snr'])
data['ra'] = data['ra']*15
paper_path = '/home/mingkang/Desktop/bootes/paper.csv'
paper = pd.read_csv(paper_path)
fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.scatter(data['ra'],data['dec'],data['deboost'])
ax.scatter(paper['R.A.'],paper['Dec.'],paper['S850'])

plt.show()