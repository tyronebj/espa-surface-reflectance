# =======================================================
#   Binary Build Layer
# =======================================================
FROM jbrinkmann/lagoon-toad:devel-1.7.0.0 as builder
LABEL maintainer="USGS EROS LSRD http://eros.usgs.gov" \
      description="ESPA scripts generating top-of-atmosphere and surface reflectance products"

WORKDIR ${SRC_DIR}
COPY . ${SRC_DIR}

# `  Install the product-formatter applications
RUN cd ${SRC_DIR}/scripts \
    && make BUILD_STATIC=yes ENABLE_THREADING=yes \
    && make install \
    && cd ${SRC_DIR}/ledaps \
    && make BUILD_STATIC=yes ENABLE_THREADING=yes \
    && make install \
    && cd ${SRC_DIR}/lasrc \
    && make BUILD_STATIC=yes ENABLE_THREADING=yes \
    && make install \
    && cd ${SRC_DIR} \
    && rm -rf *
