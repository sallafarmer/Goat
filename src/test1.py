import numpy

xs = [3, 1, 2]    # Create a list
print(xs, xs[2])  # Prints "[3, 1, 2] 2"
print(xs[-1])     # Negative indices count from the end of the list; prints "2"
xs[2] = 'foo'     # Lists can contain elements of different types
print(xs)         # Prints "[3, 1, 'foo']"
xs.append('bar')  # Add a new element to the end of the list
print(xs)         # Prints "[3, 1, 'foo', 'bar']"
x = xs.pop()      # Remove and return the last element of the list
print(x, xs)

from pprint import pprint
import rasterio
import numpy as np

path = "data/20170217_215326_0c22_1B_AnalyticMS_274501.tiff"
with rasterio.open(path) as src:
    array = src.read()

pprint(src.name)

mypath = "../data"
from os import listdir
from os.path import isfile, join
onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]

pprint(onlyfiles)

stats = []
for band in array:
    stats.append({
        'min': band.min(),
        'mean': band.mean(),
        'median': np.median(band),
        'max': band.max()})

pprint(stats)

gc = { "type": "GeometryCollection",
  "geometries": [
    { "type": "Point",
      "coordinates": [-89.33, 30.0]
    },
    { "type": "LineString",
      "coordinates": [ [-89.33, 30.30], [-89.36, 30.28] ]
    }
  ]
}

pprint(gc)

import geojson
p = geojson.Point([-92, 37])
geojs = geojson.dumps(p)
pprint(geojs)

from shapely.geometry import asShape
point = asShape(p)
pprint(point.wkt)


import shapefile
shp = shapefile.Reader("data/point/point.shp")
for feature in shp.shapeRecords():
    point = feature.shape.points[0]
rec = feature.record[0]
print(point[0], point[1], rec)


from osgeo import ogr
shp = ogr.Open("data/point/point.shp")
layer = shp.GetLayer()
for feature in layer:
    geometry = feature.GetGeometryRef()
    print(geometry.GetX(), geometry.GetY(), feature.GetField("FIRST_FLD"))

from osgeo import gdal


path = "data/20170217_215326_0c22_1B_AnalyticMS_274501.tiff"
raster = gdal.Open(path)
print(raster.RasterCount)
print(raster.RasterXSize)
print(raster.RasterYSize)

print("gdal_array")
from osgeo import gdal_array
srcArray = gdal_array.LoadFile(path)
band1 = srcArray[0]
gdal_array.SaveArray(band1, "band1.png", format="png")

print("DONE")


import operator
from osgeo import gdal, gdal_array, osr
import shapefile
try:
    import Image
    import ImageDraw
except:
    from PIL import Image, ImageDraw

# Raster image to clip
raster = path
# Polygon shapefile used to clip
shp = "data/hancock/hancock.shp"
# Name of clipped raster file(s)
output = "clip"

def imageToArray(i):
    """
    Converts a Python Imaging Library array to a gdal_array image.
    """
    a = gdal_array.numpy.fromstring(i.tobytes(), 'b')
    a.shape = i.im.size[1], i.im.size[0]
    return a

def world2Pixel(geoMatrix, x, y):
    """
    Uses a gdal geomatrix (gdal.GetGeoTransform()) to calculate
    the pixel location of a geospatial coordinate
    """
    ulX = geoMatrix[0]
    ulY = geoMatrix[3]
    xDist = geoMatrix[1]
    yDist = geoMatrix[5]
    rtnX = geoMatrix[2]
    rtnY = geoMatrix[4]
    pixel = int((x - ulX) / xDist)
    line = int((ulY - y) / abs(yDist))
    return (pixel, line)
# Load the source data as a gdal_array array
srcArray = gdal_array.LoadFile(raster)
# Also load as a gdal image to get geotransform (world file) info
srcImage = gdal.Open(raster)
geoTrans = srcImage.GetGeoTransform()
# Use pyshp to open the shapefile
r = shapefile.Reader(shp)
# Convert the layer extent to image pixel coordinates
minX, minY, maxX, maxY = r.bbox
ulX, ulY = world2Pixel(geoTrans, minX, maxY)
lrX, lrY = world2Pixel(geoTrans, maxX, minY)
# Calculate the pixel size of the new image
pxWidth = int(lrX - ulX)
pxHeight = int(lrY - ulY)
clip = srcArray[:, ulY:lrY, ulX:lrX]
# Create a new geomatrix for the image
# to contain georeferencing data
geoTrans = list(geoTrans)
geoTrans[0] = minX
geoTrans[3] = maxY
# Map points to pixels for drawing the county boundary
# on a blank 8-bit, black and white, mask image.
pixels = []
for p in r.shape(0).points:
    pixels.append(world2Pixel(geoTrans, p[0], p[1]))
rasterPoly = Image.new("L", (pxWidth, pxHeight), 1)
# Create a blank image in PIL to draw the polygon.
rasterize = ImageDraw.Draw(rasterPoly)
rasterize.polygon(pixels, 0)
# Convert the PIL image to a NumPy array
mask = imageToArray(rasterPoly)
# Clip the image using the mask
clip = gdal_array.numpy.choose(mask, (clip, 0)).astype(
                                gdal_array.numpy.uint8)
# Save ndvi as tiff
gdal_array.SaveArray(clip, "{}.tif".format(output),
                      format="GTiff", prototype=raster)
