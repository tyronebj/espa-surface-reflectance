# =======================================================
#   Binary Build Layer
# =======================================================
FROM usgseros/espa-l2qa-tools:docker-devel-3.0rc1.dev1 as builder
LABEL maintainer="USGS EROS LSRD http://eros.usgs.gov" \
      description="ESPA scripts generating top-of-atmosphere and surface reflectance products"
USER root

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
USER espadev

