# versioned rocker image: R 4.6.1 with CRAN packages from a fixed date
FROM rocker/r-ver:4.6.1

# system libraries some R packages need to install, plus pandoc for R Markdown
RUN apt-get update && apt-get install -y --no-install-recommends \
        libcurl4-openssl-dev libssl-dev libxml2-dev \
        libfontconfig1-dev libfreetype6-dev libpng-dev \
        libharfbuzz-dev libfribidi-dev libuv1-dev pandoc \
    && rm -rf /var/lib/apt/lists/*

# packages first, so this slow step is cached when only the code changes
RUN R -e "install.packages(c('KMsurv', 'survival', 'cmprsk', 'survminer', \
                             'dplyr', 'tidyr', 'ggplot2', \
                             'shiny', 'flexdashboard', 'rmarkdown', 'testthat'))"

WORKDIR /app
COPY . .

EXPOSE 3838

# default: start the dashboard; host 0.0.0.0 makes it reachable from outside the container
CMD ["R", "-e", "rmarkdown::run('dashboard.Rmd', shiny_args = list(host = '0.0.0.0', port = 3838))"]