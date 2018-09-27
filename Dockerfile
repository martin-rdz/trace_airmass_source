# Use phusion/baseimage as base image. To make your builds
# reproducible, make sure you lock down to a specific version, not
# to `latest`! See
# https://github.com/phusion/baseimage-docker/blob/master/Changelog.md
# for a list of version numbers.
FROM phusion/baseimage:latest


# Use baseimage-docker's init system.
CMD ["/sbin/my_init"]

# ...put your own build instructions here...

RUN apt-get update && add-apt-repository ppa:ubuntugis/ppa && apt-get update && apt-get install -y \
  language-pack-en openssh-server \
  build-essential make gcc zlib1g-dev git python3 python3-dev python3-pip \
  pandoc python3-setuptools\
  libatlas-base-dev gfortran \
  libjpeg8-dev libfreetype6-dev libpng12-dev \
  binutils libproj-dev libgdal-dev gdal-bin python3-numpy python3-mysqldb \
  python3-scipy

# libzmq3-dev sqlite3 libsqlite3-dev pandoc libcurl4-openssl-dev nodejs \

# set the language
ENV LANGUAGE en_US.UTF-8
ENV LANG en_US.UTF-8
ENV LC_ALL en_US.UTF-8

RUN locale-gen en_US.UTF-8
RUN dpkg-reconfigure locales

# Clean up APT when done.
RUN apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN pip3 install --upgrade pip
RUN pip3 install cython
RUN pip3 install ipython[notebook]
RUN pip3 install numpy
RUN pip3 install scipy
RUN pip3 install matplotlib
RUN pip3 install shapely
RUN pip3 install pyshp
RUN pip3 install six
RUN pip3 install affine
RUN pip3 install fastkml
RUN pip3 install toml
RUN pip3 install netcdf4
RUN pip3 install pyproj
#RUN pip3 install Mysqldb
RUN pip3 install rasterio
RUN pip3 install cartopy


#RUN rm -f /etc/service/sshd/down

# Regenerate SSH host keys. baseimage-docker does not contain any, so you
# have to do that yourself. You may also comment out this instruction; the
# init system will auto-generate one during boot.
#RUN /etc/my_init.d/00_regen_ssh_host_keys.sh

#EXPOSE 22
