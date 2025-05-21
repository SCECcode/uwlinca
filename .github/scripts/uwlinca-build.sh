#!/bin/bash

tmp=`uname -s`

if [ $tmp == 'Darwin' ]; then
##for macOS, make sure have automake/aclocal/libtoolize
  brew install automake
  brew reinstall gcc
  brew install libtool
  export PATH="/opt/homebrew/opt/libtool/libexec/gnubin:$PATH"
fi

libtoolize
aclocal -I m4
autoconf
automake --add-missing --force-missing
./configure --prefix=$UCVM_INSTALL_PATH/model/uwlinca --with-proj-libdir=$UCVM_INSTALL_PATH/lib/proj/lib --with-proj-incdir=$UCVM_INSTALL_PATH/lib/proj/include --with-tiff-libdir=$UCVM_INSTALL_PATH/lib/tiff/lib --with-tiff-incdir=$UCVM_INSTALL_PATH/lib/tiff/include --with-sqlite-libdir=$UCVM_INSTALL_PATH/lib/sqlite/lib --with-sqlite-incdir=$UCVM_INSTALL_PATH/lib/sqlite/include
make
make install

