import ipdb
import numpy as np

def parseGPX(fn):
    f = open(fn,'r')
    
    lats = []
    longs = []
    
    line = None
    while not (line == ""):
        line = f.readline()
        if line.find("<trkpt ")>-1:
            l = line.split("\"")
            lats.append(float(l[1]))
            longs.append(float(l[3]))
    
    return lats, longs

def create_waypoints(x, y):
    dists = []
    angles = []
    for i in range(len(x)-1):
        dists.append(((x[i+1]-x[i])**2 + (y[i+1]-y[i])**2)**.5)
        angles.append(np.arctan2(y[i+1]-y[i], x[i+1]-x[i]))
    cumdist = np.append([0], np.cumsum(dists))
    
    #ipdb.set_trace()
    wayx, wayy = [], []
    for dist in np.linspace(0,cumdist[-1]-0.001*cumdist[-1],100):
            ind = np.nonzero(cumdist>dist)[0][0]-1
            wayx.append(np.cos(angles[ind])*(dist-cumdist[ind]) + x[ind])
            wayy.append(np.sin(angles[ind])*(dist-cumdist[ind]) + y[ind])
    
    return wayx, wayy

def get_waypoints_from_gpx(fn):
    lats, lons = parseGPX(fn)
    return create_waypoints(lats, lons)