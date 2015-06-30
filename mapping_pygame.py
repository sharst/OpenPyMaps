import numpy as np
import glob, ipdb, math, os
import urllib
import random
from scipy.ndimage.filters import gaussian_filter
from scipy.ndimage import zoom as image_zoom
import CoordinateConverter
from tiledownloader import Downloader
import pygame
from collections import deque

# XXX: The conversions into pixels should really be done within the map object

class Map():
    def __init__(self, center_coo, size, zoom):
        self.downloader = Downloader(self.tilecallback)
        self.size = size
        self.display = pygame.display.set_mode(self.size)
        # Calculate how many tiles we need to cache.
        max_tiles = (int(self.size[0]/256.)+1) * (int(self.size[1]/256.)+1) 
        self.tile_cache = deque([], maxlen = max_tiles)
        if isinstance(center_coo, list):
            center_coo = Coordinate().from_lat_long(center_coo[0], center_coo[1])
        elif isinstance(center_coo, Coordinate):
            pass
        else:
            raise ValueError("Please provide a Coordinate instance or a lat/long pair")
        
        self.set_map(center_coo, zoom)
        
        # Register event handler
        self.cid = self.fig.canvas.mpl_connect('button_press_event', self.onclick)
        
        self.overlays = {}
        self.layers = {}
        self.draw()
        
    def tilecallback(self, tile):
        (x,y) = self.deg2pix(*num2deg(tile.x, tile.y, tile.zoom))
        pass
        #XXX Do something!
    
    def onclick(self, event):
        print 'button=%d, x=%d, y=%d, xdata=%f, ydata=%f'%(
            event.button, event.x, event.y, event.xdata, event.ydata)
        lat, lon = num2deg(event.xdata/256.+self.x_num, event.ydata/256.+self.y_num, self.zoom)
        self.set_map(Coordinate().from_lat_long(lat, lon), self.zoom)
        self.draw()
    
    def set_map(self, center_coo, zoom):
        # XXX: If there is a heatmap on top, it needs to be moved, too.
        self.zoom = zoom
        self.center_coo = center_coo
        (x, y, zoom) = deg2num(center_coo.get_lat_long()[0], center_coo.get_lat_long()[1], zoom)
                
        x_l, y_l = x-(self.size[0]/2)/256., y-(self.size[1]/2)/256.
        x_r, y_r = x+(self.size[0]/2)/256., y+(self.size[1]/2)/256.
        
        for x in range(int(x_l), int(x_r)+1):
            for y in range(int(y_l), int(y_r)+1):
                self.downloader.request_download(x, y, zoom)

        self.x_num = x_l
        self.y_num = y_l
    
    def draw(self):
        self.ax.clear()
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
        self.fig.canvas.draw()
    
    
    def add_point_cloud(self, name, pointcloud):
        self.overlays[name] = pointcloud
        self.draw()
    
    def remove_point_cloud(self, name):
        try:
            self.overlays.pop(name)
            self.draw()
        except KeyError:
            print "No Pointcloud found with name %s " %name
    
    def add_heatmap(self, name, heatmap):
        self.layers[name] = heatmap
        self.draw()
        
    def remove_heatmap(self, name):
        self.layers.pop(name)
        self.draw()
    
    """
    Add a marker to be drawn on this map.
    @param name: a unique name for this marker
    @param marker: A mapping.Marker instance
    """
    def add_marker(self, name, marker):
        self.overlays[name] = marker
        self.draw()

    """
    Remove the marker of that name from the map.
    """
    def remove_marker(self, name):
        self.overlays.pop(name)
        self.draw()
    
    """ Convert latitude, longitude coordinates into pixels for this map 
        @param lat: Latitude in degreees
        @param lon: Longitude in degrees 
    """
    def deg2pix(self, lat, lon):
        # Convert endpoint to appropriate map projection
        (xtile2, ytile2, zoom) = deg2num(lat, lon, self.zoom)
        # Convert endpoint to figure pixels           
        x = (xtile2-self.x_num)*256
        y = (ytile2-self.y_num)*256
        return (x,y)
    
    def set_zoom(self, zoom):
        self.set_map(self.center_coo, zoom)
        self.draw()
    
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
    def __init__(self, map, xdata, ydata, color = "r", marker = "o", linestyle = "-"):
        self.lats = xdata
        self.longs = ydata
        self.color = color
        self.marker = marker
        self.ls = linestyle
    
    def get_patches(self, map):
        x = [map.deg2pix(lat, lon)[0] for (lat, lon) in zip(self.lats, self.longs)]
        y = [map.deg2pix(lat, lon)[1] for (lat, lon) in zip(self.lats, self.longs)]
        self.patch = mpl.lines.Line2D(x, y, color = self.color, linestyle = self.ls, marker = self.marker)
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
        marker = Marker([60.291983,-44.280996], 20, "Home", 20, [10,30])
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
            