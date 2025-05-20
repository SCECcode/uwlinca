#!/usr/bin/env python

##
#  Builds the data files in the expected format from 
#
# from >> easting(Km) northing(Km) depth(Km) vp(km/s)
#
# to >>  /** P-wave velocity in km/s per second */
#        double vp;
# depth is in increment of 1000m,
#
#The columns of the file are: x (km) y (km) z (km) vp (km/s), where 
#x, y, and z are utm coordinates in km and the columns increase in x.y.z order. 
#Grid spacing is 1 km in each spatial dimension. 

import getopt
import sys
import subprocess
import struct
import array
import ssl
import certifi
from urllib.request import urlopen

if sys.version_info.major >= (3) :
  from urllib.request import urlopen
else:
  from urllib2 import urlopen

## at UWLINCA/MOD.finer1

model = "UWLINCA"

dimension_x = 0
dimension_y = 0 
dimension_z = 0

def usage():
    print("\n./make_data_files.py\n\n")
    sys.exit(0)

def download_urlfile(url,fname):
  print("\ndata file:",url,"\n")
  try:
    response = urlopen(url)
    CHUNK = 16 * 1024
    with open(fname, 'wb') as f:
      while True:
        chunk = response.read(CHUNK)
        if not chunk:
          break
        f.write(chunk)
  except:
    e = sys.exc_info()[0]
    print("Exception retrieving and saving model datafiles:",e)
    raise
  return True

def download_urlfile2(url, fname):
    print("\ndata file:", url, "\n")
    try:
        context = ssl.create_default_context(cafile=certifi.where())
        response = urlopen(url, context=context)
        CHUNK = 16 * 1024
        with open(fname, 'wb') as f:
            while True:
                chunk = response.read(CHUNK)
                if not chunk:
                  break
                f.write(chunk)
    except Exception as e:
        print("Exception retrieving and saving model datafiles:", e)
        raise
    return True


def main():

    # Set our variable defaults.
    path = ""
    mdir = ""

    try:
        fp = open('./config','r')
    except:
        print("ERROR: failed to open config file")
        sys.exit(1)

    ## look for model_data_path and other varaibles
    lines = fp.readlines()
    for line in lines :
        if line[0] == '#' :
          continue
        parts = line.split('=')
        if len(parts) < 2 :
          continue;
        variable=parts[0].strip()
        val=parts[1].strip()

        if (variable == 'model_data_path') :
            path = val + '/' + model
            continue
        if (variable == 'model_dir') :
            mdir = "./"+val
            continue
        if (variable == 'nx') :
            dimension_x = int(val)
            continue
        if (variable == 'ny') :
            dimension_y = int(val)
            continue
        if (variable == 'nz') :
            dimension_z = int(val)
            continue
        continue
    if path == "" :
        print("ERROR: failed to find variables from config file")
        sys.exit(1)

    fp.close()

    print("\nDownloading model file\n")

    fname="./"+"MOD.finer1"
    url = path + "/" + fname
    download_urlfile2(url,fname)

    subprocess.check_call(["mkdir", "-p", mdir])

    # Now we need to go through the data files and put them in the correct
    # format for LSU_IV. More specifically, we need a Vp.dat

    fvp = open("./MOD.finer1");
    f_vp = open("./uwlinca/vp.dat", "wb")

    vp_arr = array.array('f', (-1.0,) * (dimension_x * dimension_y * dimension_z))

    print ("dimension is", (dimension_x * dimension_y * dimension_z))

    vp_nan_cnt = 0
    vp_total_cnt =0;

    x_pos=0;
    y_pos=0;
    z_pos=0;

    for line in fvp:
        arr = line.split()

        vp = -1.0
        in_x = float(arr[0])
        in_y = float(arr[1])
        in_z = float(arr[2])
        tmp = arr[3]

        if( tmp != "NaN" ) :
           vp = float(tmp)
           vp = vp * 1000.0;
        else:
           vp_nan_cnt = vp_nan_cnt + 1

        vp_total_cnt = vp_total_cnt + 1

        loc =z_pos * (dimension_y * dimension_x) + (y_pos * dimension_x) + x_pos
        vp_arr[loc] = vp

        x_pos = x_pos + 1
        if(x_pos == dimension_x) :
          x_pos = 0;
          y_pos = y_pos+1
          if(y_pos == dimension_y) :
            y_pos=0;
            z_pos = z_pos+1
            if(z_pos == dimension_z) :
              print ("All DONE")

    vp_arr.tofile(f_vp)

    fvp.close()
    f_vp.close()

    print("Done! with NaN(", vp_nan_cnt, ") total(", vp_total_cnt,")")


if __name__ == "__main__":
    main()

