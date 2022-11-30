
#FROM ubuntu:latest
FROM ubuntu:22.04


ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get update && apt-get install -y \
  language-pack-en openssh-server vim software-properties-common \
  build-essential make gcc g++ zlib1g-dev git python3 python3-dev python3-pip \
  pandoc python3-setuptools imagemagick\
  gfortran autoconf libtool automake flex bison cmake git-core \
  libjpeg8-dev libfreetype6-dev libhdf5-serial-dev \
  libeccodes0 libeccodes-data libeccodes-dev \
  libnetcdf-c++4 libnetcdf-c++4-dev libnetcdff-dev \
  binutils  python3-numpy python3-mysqldb \
  python3-scipy python3-sphinx libedit-dev unzip curl wget \
  libatlas-base-dev libpng-dev \
  libnetcdff-dev libopenjp2-7-dev gfortran make unzip git cmake wget \
  libproj-dev libgdal-dev gdal-bin 
  
  
# replaced 'libpng12-dev' by libpng-dev
#RUN add-apt-repository ppa:ubuntugis/ppa \
#  && apt-get update \
#  && apt-get install -y libatlas-base-dev libpng-dev \
#     libproj-dev libgdal-dev gdal-bin  

#RUN add-apt-repository 'deb http://security.ubuntu.com/ubuntu xenial-security main' \
#  && apt-get update \
#  && apt-get install -y libjasper1 libjasper-dev libeccodes-tools libeccodes-dev

#ENV HTTP https://confluence.ecmwf.int/download/attachments/45757960
#ENV ECCODES eccodes-2.9.2-Source


RUN mkdir source_builds && cd source_builds && mkdir eccodes && cd eccodes \
  && wget https://confluence.ecmwf.int/download/attachments/45757960/eccodes-2.27.0-Source.tar.gz?api=v2 \
  && tar -xzf eccodes-2.27.0-Source.tar.gz?api=v2 \
  && mkdir build && cd build \
  && cmake -DCMAKE_INSTALL_PREFIX=/usr/ ../eccodes-2.27.0-Source \
  #&& make && ctest --output-on-failure && make install
  && make && make install

#RUN cd /tmp && wget --output-document=$ECCODES.tar.gz $HTTP/$ECCODES.tar.gz?api=v2 \
#  && tar -xzf $ECCODES.tar.gz \
#  && mkdir -p /tmp/$ECCODES/build \
#  && cd /tmp/$ECCODES/build \
#  && cmake \
#      -DCMAKE_INSTALL_PREFIX=/usr \
#      -DCMAKE_BUILD_TYPE=production \
#      -DENABLE_PYTHON=1 -DENABLE_MEMFS=1 \
#      -DENABLE_ECCODES_THREADS=1 .. \
#  && make -j2 \
#  && make install \
#  && cd / \
#  && rm -rf /tmp/$ECCODES*

#
# Download, modify and compile flexpart 10
#
#RUN mkdir flex_src && cd flex_src \
#  && wget https://www.flexpart.eu/downloads/66 \
#  && tar -xvf 66 \
#  && rm 66 \
#  && cd flexpart_v10.4_3d7eebf/src \
#  && cp makefile makefile_local \
#  && sed -i '74 a INCPATH1 = /usr/include\nINCPATH2 = /usr/include\nLIBPATH1 = /usr/lib\n F90 = gfortran' makefile_local \
#  && sed -i 's/LIBS = -lgrib_api_f90 -lgrib_api -lm -ljasper $(NCOPT)/LIBS = -leccodes_f90 -leccodes -lm -ljasper $(NCOPT)/' makefile_local \
#  && make -f makefile_local

RUN mkdir flex_src && cd /flex_src \
  && git clone https://www.flexpart.eu/gitmob/flexpart  --branch dev --single-branch \
  && cd flexpart/src \
  && cp makefile makefile_local \
  && sed -i '74 a INCPATH1 = /usr/include\nINCPATH2 = /usr/include\nLIBPATH1 = /usr/lib\n F90 = gfortran' makefile_local \
#  && sed -i 's/LIBS = -lgrib_api_f90 -lgrib_api -lm -ljasper $(NCOPT)/LIBS = -leccodes_f90 -leccodes -lm -ljasper $(NCOPT)/' makefile_local
  && sed -i 's/LIBS = -lgrib_api_f90 -lgrib_api -lm -ljasper $(NCOPT)/LIBS = -leccodes_f90 -leccodes -lm -ljasper $(NCOPT)/' makefile_local \
  && make -f makefile_local

#ENV PATH /flex_src/flexpart_v10.4_3d7eebf/src/:$PATH
ENV PATH /flex_src/flexpart/src/:$PATH

# # set the language
# ENV LANGUAGE en_US.UTF-8
# ENV LANG en_US.UTF-8
# ENV LC_ALL en_US.UTF-8

# RUN locale-gen en_US.UTF-8
# RUN dpkg-reconfigure locales

RUN pip3 install --upgrade pip
RUN pip3 install cython
RUN pip3 install ipython[notebook]
RUN pip3 install numpy
RUN pip3 install scipy
RUN pip3 install numba
RUN pip3 install matplotlib==2.0.2
RUN pip3 install shapely --no-binary shapely
RUN pip3 install pyshp
RUN pip3 install six
RUN pip3 install eccodes
RUN pip3 install xarray
RUN pip3 install cfgrib
RUN pip3 install affine
RUN pip3 install fastkml
RUN pip3 install toml
RUN pip3 install netcdf4
RUN pip3 install pyproj
#RUN pip3 install Mysqldb
RUN pip3 install rasterio
RUN pip3 install cartopy
RUN pip3 install --upgrade Sphinx
RUN pip3 install sphinx_rtd_theme
RUN pip3 install recommonmark
RUN pip3 install bcolz-zipline
RUN pip3 install pyhdf

EXPOSE 8890
