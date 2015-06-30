import threading
import random
import os
import urllib
import pygame
import time


class MapTile(object):
    def __init__(self, x, y, zoom):
        self.x = x
        self.y = y
        self.zoom = zoom
        self.data = []
        
    
class Downloader(threading.Thread):
    def __init__(self, callback):
        super(Downloader, self).__init__()
        self.callback = callback
        self.download_list = []
        # XXX This should of course be a temp folder at some point
        if not os.path.exists("tiles/"):
            os.mkdir("tiles/")
        self.running = True
        
    def request_download(self,x,y,zoom):
        self.download_list.append(MapTile(x,y,zoom))
        
    def run(self):
        while self.running:
            if len(self.download_list) > 0:
                tile = self.download_list.pop(0)
                # Gets a tile from cache or if not present downloads from osm server
                tilename = str(tile.zoom) + '/' + str(tile.x) + '/' +str(tile.y) + '.png'
                
                if (os.path.isfile('tiles/'+tilename.replace('/','_'))):
                    print "Fetching image from file: " + tilename.replace('/','_')
                else:
                    URL = 'http://' + random.choice(['a','b','c'])+ '.tile.openstreetmap.org/' + tilename
                    print "Fetching image from URL: " + URL
                    urllib.urlretrieve(URL, filename = 'tiles/'+tilename.replace('/','_'))
                
                tile.data = pygame.image.load('tiles/'+tilename.replace('/','_'))   
                self.callback(tile)
            else:
                time.sleep(.1)
            
    def shutdown(self):
        self.running = False