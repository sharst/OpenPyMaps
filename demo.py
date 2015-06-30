import mapping
import gpxparsing
import pylab as plt

mm = mapping.Map([52.28764005, 7.984405159], [600,600], 16)
x, y = gpxparsing.get_waypoints_from_gpx("SeeRoute.gpx")
mm.add_point_cloud("Track", mapping.Pointcloud(mm, x, y, "r", "o"))
raw_input()