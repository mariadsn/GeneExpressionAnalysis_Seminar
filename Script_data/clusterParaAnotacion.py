# crear ficheros para clusters 


import os

path = os.path.dirname(os.path.realpath('GenclusterParaAnotacion.pyr'))


lClusterG = {}
clusterD = open(path +'/TumorBreast_DEGs_table_GLay.csv','r')
for line in clusterD:
    line = line.strip()
    if not line.startswith('"__glayCluster"'):
        spl = line.split(',')
        clstr = spl[0].replace('"','')
        name = spl[1].replace('"','')
        if not clstr in lClusterG:
            lClusterG[clstr] = []
        if not name in lClusterG[clstr]:
            lClusterG[clstr].append(name)
clusterD.close()
Good_clust = []
for n in lClusterG:
    if len(lClusterG[n]) >10:
        Good_clust.append(n)
print(len(Good_clust))

for elem in lClusterG:
    if len(lClusterG[elem])>10:
        f = open(path + '/Cluster_' + str(elem) + '.txt','w')
        for l in lClusterG[elem]:
            f.write(l + '\n')
        f.close()

