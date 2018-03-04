from pprint import pprint
import rasterio
import numpy as np
import numpy.ma as ma
import math
import matplotlib.pyplot as plt


def getMeta(path):
    with rasterio.open(path) as src:
        meta = src.meta
    return meta

def readTiff(path):
    with rasterio.open(path) as src:
        profile = src.profile
        pprint(src.shape)
        print(src.crs)
        print(src.count)
        print(src.indexes)
        print(src.bounds)
        print(src.transform)
        array = src.read()
    pprint(src.name)
    return array

def genZoneMasks(zones):
    zone = np.array(zones[0])
    cleaned = np.ma.array(zone, mask=np.isnan(zone))
    uniq = np.unique(cleaned)
    zoneArray = []
    i = 1
    for y in uniq:
        m1 = ma.masked_not_equal(zone, y)
        np.putmask(m1, m1 == y, 0.0)
        zoneArray.append(m1)
        i = i + 1
    k = 0
    t = 0
    for x in np.nditer(zone):
        t = t+1
        if not math.isnan(x) :
            k = k+1
    print k, t
    print uniq
    return zoneArray

def genMaskedBand(band, zoneArray):
    bandArr = np.array(band)
    maskedBandArr = []
    for x in zoneArray:
        maskBand = np.ma.array(bandArr, mask=x)
        maskedBandArr.append(maskBand)
    return maskedBandArr


def calcNdvi(red, nir):
    """Calculate NDVI."""
    return ((nir - red) / (nir + red))*10000

def generateNDVIPerZone(tiffFile, zoneMapFile):
    ndviZones = dict()
    ndviCalcs = dict()
    zones = readTiff(zoneMapFile)
    original = readTiff(tiffFile)
    zoneArray = genZoneMasks(zones)

    numZones = len(zoneArray)
    """BGRN"""
    maskBandArr0 = genMaskedBand(original[0], zoneArray)
    maskBandArr1 = genMaskedBand(original[1], zoneArray)
    maskBandArr2 = genMaskedBand(original[2], zoneArray)
    maskBandArr3 = genMaskedBand(original[3], zoneArray)

    zoneMapFileName = zoneMapFile.split('/')[3]
    tiffFileName = tiffFile.split('/')[2] + zoneMapFileName

    i = 0
    for x in range(len(zoneArray)-1):
        ndvi = calcNdvi(maskBandArr2[i], maskBandArr3[i])
        ndviZones[x] = ndvi
        ndviCalcs[x] = (np.mean(ndvi), np.std(ndvi))
        meta = getMeta(tiffFile)
        meta.update(dtype=rasterio.uint16, count=1)
        with rasterio.open("../data/output/{}_ndvi{}.tiff".format(tiffFileName, i), 'w', **meta) as out:
            out.write_band(1, ndvi)
        plt.imsave("../data/output/{}_ndvi{}.png".format(tiffFileName, i), ndvi, cmap=plt.cm.summer)
        i = i+1

    meta.update(dtype=rasterio.uint16, count=4)
    with rasterio.open("../data/output/{}_orig.tiff".format(tiffFileName), 'w', **meta) as out:
        out.write_band(1, original[0])
        out.write_band(2, original[1])
        out.write_band(3, original[2])
        out.write_band(4, original[3])

    plt.imsave("../data/output/{}_originalB.png".format(tiffFileName), original[0], cmap=plt.cm.summer)
    plt.imsave("../data/output/{}_originalG.png".format(tiffFileName), original[1], cmap=plt.cm.summer)
    plt.imsave("../data/output/{}_originalR.png".format(tiffFileName), original[2], cmap=plt.cm.summer)
    plt.imsave("../data/output/{}_originalN.png".format(tiffFileName), original[3], cmap=plt.cm.summer)

    for k, v in ndviCalcs.items():
        print k
        pprint(v)

    # for k,v in ndviZones.items():
    #    print k
    #    pprint(v)

    return (ndviCalcs, ndviZones)


zoneMapFile = "../data/(controller)_2017-08-30-2017-09-05__274501__Harvest/2017-08-30-2017-09-05__274502__Harvest__ZoneMap.tiff"
tiffFile = "../data/20170217_215326_0c22_1B_AnalyticMS_274501.tiff"

#(ndviCalcs, ndviZones) = generateNDVIPerZone(tiffFile, zoneMapFile)

zoneMapFile = "../data/(controller)_2017-09-13-2017-09-28__297113__Harvest/2017-09-13-2017-09-28__297114__Harvest__ZoneMap.tiff"
zoneMapFile2 = "../data/(controller)_2017-09-13-2017-09-28__297113__Harvest/2017-09-13-2017-09-28__297115__Harvest__ZoneMap.tiff"
zoneMapFile3 = "../data/(controller)_2017-09-13-2017-09-28__297113__Harvest/2017-09-13-2017-09-28__428989__Harvest__ZoneMap.tiff"
tiffFile = "../data/20170223_165518_0e30_1B_AnalyticMS_297113.tiff"

print zoneMapFile
(ndviCalcs, ndviZones) = generateNDVIPerZone(tiffFile, zoneMapFile)

print zoneMapFile2

(ndviCalcs2, ndviZones) = generateNDVIPerZone(tiffFile, zoneMapFile2)

print zoneMapFile3
(ndviCalcs3, ndviZones) = generateNDVIPerZone(tiffFile, zoneMapFile3)


pprint("Completed")



