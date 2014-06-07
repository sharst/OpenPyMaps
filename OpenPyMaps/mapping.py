import pylab as pp
import numpy as np
import glob, ipdb, math, os
import urllib
import random
import matplotlib.image as mpimg
import matplotlib as mpl
from scipy.ndimage.filters import gaussian_filter
from scipy.ndimage import zoom as image_zoom
import CoordinateConverter

# XXX: The conversions into pixels should really be done within the map object

class Map():
    def __init__(self, center_coo, size, zoom):
        if isinstance(center_coo, list):
            center_coo = Coordinate().from_lat_long(center_coo[0], center_coo[1])
        elif isinstance(center_coo, Coordinate):
            pass
        else:
            raise ValueError("Please provide a Coordinate instance or a lat/long pair")
        
        # XXX This should of course be a temp folder at some point
        if not os.path.exists("tiles/"):
            os.mkdir("tiles/")
            
        (x, y, zoom) = deg2num(center_coo.get_lat_long()[0], center_coo.get_lat_long()[1], zoom)
        
        x, y = int(x), int(y)
        #asdfsdfadfasf
        data = None
        for x in range(x, x+1):
            datacolumn = None
            for y in range(y_l, y_r+1):
                tmp = getTile((x,y,zoom))
                if datacolumn is None:
                    datacolumn = tmp
                else:
                    datacolumn = np.append(datacolumn,tmp, axis = 0)
            if data is None:
                data = datacolumn
            else:
                data = np.append(data, datacolumn, axis = 1)
        
        # The map, tilenumbers of left upper and lower right tile
        self.ax = None
        self.map = data
        self.xtilerange=(x_l,x_r)
        self.ytilerange=(y_l,y_r)
        self.zoom = zoom
        
        left_ext = num2deg(x_l, y_l, zoom)
        right_ext = num2deg(x_r+1, y_r+1, zoom)
        self.extent = [left_ext[0], left_ext[1],right_ext[0],right_ext[1]]
        
        self.overlays = {}
        self.layers = {}
    
    def draw(self):
        if self.ax is None:
            self.fig = pp.figure("Map")
            #mng = pp.get_current_fig_manager()
            #mng.resize(*mng.window.maxsize())
            self.ax = pp.subplot(111)
            self.ax.imshow(self.map)
            
            markers = [elem for marker in self.overlays.values() 
                       for elem in marker.get_patches(self)]
            
            for marker in markers:
                self.ax.add_artist(marker)
            
            for hm in self.layers.values():
                self.ax.imshow(hm.get_heatmap(self), cmap = hm.cmap, alpha = hm.opacity)
            
            self.ax.set_yticks([])
            self.ax.set_xticks([])
            
            self.ax.set_xlim((0, self.map.shape[1]))
            self.ax.set_ylim((self.map.shape[0], 0))
        else:
            self.fig.canvas.draw()
    
    
    def add_point_cloud(self, name, pointcloud):
        self.overlays[name] = pointcloud
    
    def remove_point_cloud(self, name):
        self.overlays.pop(name)
    
    def add_heatmap(self, name, heatmap):
        self.layers[name] = heatmap
        
    def remove_heatmap(self, name):
        self.layers.pop(name)
    
    """
    Add a marker to be drawn on this map.
    @param name: a unique name for this marker
    @param marker: A mapping.Marker instance
    """
    def add_marker(self, name, marker):
        self.overlays[name] = marker

    """
    Remove the marker of that name from the map.
    """
    def remove_marker(self, name):
        self.overlays.pop(name)
    
    """ Convert latitude, longitude coordinates into pixels for this map 
        @param lat: Latitude in degreees
        @param lon: Longitude in degrees 
    """
    def deg2pix(self, lat, lon):
        # Convert endpoint to appropriate map projection
        (xtile2, ytile2, zoom) = deg2num(lat, lon, self.zoom)
        # Convert endpoint to figure pixels           
        x = (xtile2-self.xtilerange[0])*256
        y = (ytile2-self.ytilerange[0])*256
        return (x,y)
    
""" Gives the tile number of osm-tiles for lat and long + zoom """
def deg2num(lat_deg, lon_deg, zoom):
    lat_rad = math.radians(lat_deg)
    n = 2.0 ** zoom
    xtile = (lon_deg + 180.0) / 360.0 * n
    ytile = (1.0 - math.log(math.tan(lat_rad) + (1 / math.cos(lat_rad))) / math.pi) / 2.0 * n
    return (xtile, ytile, zoom)

""" Gives the coordinates of upper left corner of a tile. For lower right use xtile+1, ytile+1 """
def num2deg(xtile, ytile, zoom):
    n = 2.0 ** zoom
    lon_deg = xtile / n * 360.0 - 180.0
    lat_rad = math.atan(math.sinh(math.pi * (1 - 2 * ytile / n)))
    lat_deg = math.degrees(lat_rad)
    return (lat_deg, lon_deg)

def get_zone_number(self, lat, lon):
    converter = CoordinateConverter.Converter([lat, long], type=CoordinateConverter.LATLON)
    return converter.zone_letter, converter.zone_number

# Gets a tile from cache or if not present downloads from osm server
def getTile((xtile,ytile,zoom)):
    tilename = str(zoom) + '/' + str(xtile) + '/' +str(ytile) + '.png'
    
    URL = 'http://' + random.choice(['a','b','c'])+ '.tile.openstreetmap.org/' + tilename
    if (os.path.isfile('tiles/'+tilename.replace('/','_'))):
        print "Fetching image from file: " + tilename.replace('/','_')
    else:
        print "Fetching image from URL: " + URL
        urllib.urlretrieve(URL, filename = 'tiles/'+tilename.replace('/','_'))
        
    return mpimg.imread('tiles/'+tilename.replace('/','_'))


class Coordinate():
    def __init__(self):
        self.x = 0.
        self.y = 0.
        self.converter = None 
        
    def from_gauss_krueger(self, right, height):
        if self.converter is None:
            self.converter = CoordinateConverter.Converter([right, height], type=CoordinateConverter.GAUSSKRUGER)
        
        x,y = self.converter.GK2UTM(right, height)
        self.x = x
        self.y = y
        return self
    
    def from_lat_long(self, lat, long):
        if self.converter is None:
            self.converter = CoordinateConverter.Converter([lat, long], type=CoordinateConverter.LATLON)
        x,y = self.converter.LatLong2UTM(lat, long)
        self.x = x
        self.y = y
        return self
    
    def from_UTM(self, x, y, zone_letter, zone_number):
        if self.converter is None:
            self.converter = CoordinateConverter.Converter()
            self.converter.zone_letter = zone_letter
            self.converter.zone_number = zone_number
        self.x = x
        self.y = y
        return self
    
    def get_lat_long(self):
        return self.converter.UTM2LatLong(self.x, self.y)
    
    def get_zone_letter_number(self):
        return self.converter.zone_letter, self.converter.zone_number
    
    def get_UTM(self):
        return self.x, self.y 


class Heatmap():
    """
    Create a heatmap from some measured data. Such a heatmap can be added as an overlay to a map.
    
    @param lats: The latitudes at which the data was measured
    @param longs: The longitudes at which the data was measured
    @param data: The numerical data that was measured 
    @param cmap: The colormap to be used when rendering this heatmap
    """
    def __init__(self, lats, longs, data, cmap = pp.cm.jet, opacity = .6):
        self.lats = np.array(lats)
        self.longs = np.array(longs)
        self.data = np.array(data)
        self.cmap = cmap
        self.opacity = opacity
    
    """
    Get an overlay representation that fits on the provided map
    """
    def get_heatmap(self, map):
        [x1,y1,x2,y2] = map.extent
        (num_pix_y, num_pix_x, _) = map.map.shape
        bins_x = np.linspace(x1,x2,num_pix_x/10.)
        bins_y = np.linspace(y1,y2,num_pix_y/10.)
        
        H = np.zeros((len(bins_y),len(bins_x)))
        
        for i in range(len(bins_x)-1):
            for j in range(len(bins_y)-1):
                H[j,i] = np.sum(self.data[self.values_between(self.longs,bins_y[j],bins_y[j+1]) & 
                                          self.values_between(self.lats,bins_x[i],bins_x[i+1])])
        
        H_filt = gaussian_filter(H,3)
        H = image_zoom(H_filt, 10.)
        
        return H        
    
    def values_between(self, x, x1, x2):
        if x1<x2:
            return (x1<x) & (x<=x2)
        elif x2<x1:
            return (x2<x) & (x<=x1)
        else:
            return False

class Pointcloud():
    # XXX: Data should be passed as coordinates!
    def __init__(self, map, xdata, ydata, color = "r", marker = "o"):
        self.lats = xdata
        self.longs = ydata
        self.color = color
        self.marker = marker
        x = [map.deg2pix(lat, lon)[0] for (lat, lon) in zip(self.lats, self.longs)]
        y = [map.deg2pix(lat, lon)[1] for (lat, lon) in zip(self.lats, self.longs)]
        self.patch = mpl.lines.Line2D(x,y, color = self.color, marker = self.marker)
    
    def get_patches(self, map):
        return [self.patch]

    def change_data(self, map, xdata, ydata):
        self.lats = xdata
        self.longs = ydata
        x = [map.deg2pix(lat, lon)[0] for (lat, lon) in zip(self.lats, self.longs)]
        y = [map.deg2pix(lat, lon)[1] for (lat, lon) in zip(self.lats, self.longs)]
        self.patch.set_data(x, y)

        
    


class Marker():
    # XXX: Muss noch so umgeschrieben werden, dass man die position veraendern kann. (Get patches muss referenzen zurueckgeben!)
    """
    Create a new marker. At the moment, a marker is always circular, with a given size. If you want a marker
    to represent a direction, you can additionally specify the direction parameter, which adds an arrow to the marker.
    If the certainty parameter is given, a circle segment is added to represent the certainty of the direction.
    
    @param coo: The coordinates of the marker touple or list of (lat, lon), in degrees
    @param size: The size of the marker, in pixels
    @param text: Optionally, text to display along with this marker 
    @param direction: Optionally, direction to display by the means of an arrow (in degrees, clockwise from north)
    @param certainty: Optionally, certainty of the direction to display along the arrow (touple of absolute angles clockwise from north)
    @param length: length of the arrow
    @param color: color of arrow and marker   
    
    Example:
        marker = Marker("Home", mapping.Marker([60.291983,-44.280996], 20, "Home", 20, [10,30])
    """
    def __init__(self, coo, size, text=None, direction = None, certainty = None, length = 100, color = 'r'):
        self.lat = coo[0]
        self.lon = coo[1]
        self.size = size
        self.length = length
        self.color = color
        self.direction = direction
        self.certainty = certainty
        self.text = text

        if direction is not None:
            # Calculate endpoint of an arrow pointing at given angle
            start = np.radians(coo)
            lat2 = np.arcsin(np.sin(start[0])*np.cos(0.1/6371.01) + 
                             np.cos(start[0])*np.sin(0.1/6371.01)*np.cos(np.radians(direction)));
            lon2 = start[1] + np.arctan2(np.sin(np.radians(direction))*np.sin(0.1/6371.01)*np.cos(start[0]), 
                                       np.cos(0.1/6371.01)-np.sin(start[0])*np.sin(lat2));
            
            self.lat2, self.lon2 = lat2, lon2

        
    def get_patches(self, map):
        patches  = []
        # Figure out pixel values for given coordinates
        (x,y) = map.deg2pix(self.lat, self.lon)
        patches.append(pp.Circle((x,y),self.size,color=self.color))
        
        if self.direction is not None:
            x2,y2 = map.deg2pix(np.degrees(self.lat2), np.degrees(self.lon2))
            arrow = pp.arrow(x, y, x2-x, y2-y, shape='full', lw = 3, length_includes_head = True, head_width = 20, color = self.color)
            patches.append(arrow)
        
        if self.certainty is not None:
            wedge = mpl.patches.Wedge((x, y), 60, self.certainty[0]-90, self.certainty[1]-90, alpha = .7, color = 'b')
            patches.append(wedge)
        
        if self.text is not None:
            text = pp.text(x, y, self.text)
            patches.append(text)
        
        return patches
            