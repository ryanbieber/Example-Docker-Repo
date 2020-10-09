FROM rocker/tidyverse:4.0.0

# update some packages
RUN apt-get update \
  && apt-get install -y --no-install-recommends \
    file \
    libxml2-dev \
    libpq-dev \
    zlib1g-dev \
    default-jdk \
    default-jre \
    libbz2-dev \
    libpcre3-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev


# copy the setup script, run it, then delete it
COPY src/setup.R /
RUN Rscript setup.R && rm setup.R

RUN mkdir /src

# copy all the other R files.
COPY src /src

# copy all the data over if you have any hard coded data
COPY data /src/data


WORKDIR /src


## run script
ENTRYPOINT Rscript -e "source('./main.R')"

