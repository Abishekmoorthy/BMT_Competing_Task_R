# BMT DATASET: PREPROCESSING
# Run this file first to load and prepare bmt_clean for analysis.

required_packages <- c("KMsurv", "tidyverse")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}


library(KMsurv) # This library is for lodaing the inbuilt R dataset
library(tidyverse) 
data("bmt", package = "KMsurv")
bmt_original <- bmt
bmt_clean <- bmt_original


dim(bmt_clean) # printing no of rows and columns
names(bmt_clean) # column names
head(bmt_clean) # first 6 rows
tail(bmt_clean) # last 6 rows

glimpse(bmt_clean) # dataset structure

summary(bmt_clean) # summary statistics
sum(duplicated(bmt_clean))
sum(is.na(bmt_clean))

bmt_clean <- bmt %>%
  mutate(
    group = factor(
      group,
      levels = c(1, 2, 3),
      labels = c("ALL", "AML low risk", "AML high risk")
    )
  )

bmt_clean <- bmt_clean %>%
  mutate(
    z3 = factor(
      z3,
      levels = c(0, 1),
      labels = c("Female", "Male")
    )
  )

bmt_clean <- bmt_clean %>%
  mutate(
    z4 = factor(
      z4,
      levels = c(0, 1),
      labels = c("Female", "Male")
    )
  )

bmt_clean <- bmt_clean %>%
  mutate(
    z5 = factor(
      z5,
      levels = c(0, 1),
      labels = c("Negative", "Positive")
    )
  )

bmt_clean <- bmt_clean %>%
  mutate(
    z6 = factor(
      z6,
      levels = c(0, 1),
      labels = c("Negative", "Positive")
    )
  )

bmt_clean <- bmt_clean %>%
  mutate(
    z8 = factor(
      z8,
      levels = c(0, 1),
      labels = c("Other", "FAB grade 4 or 5 and AML")
    )
  )

bmt_clean <- bmt_clean %>%
  mutate(
    z9 = factor(z9)
  )

bmt_clean <- bmt_clean %>%
  mutate(
    z10 = factor(
      z10,
      levels = c(0, 1),
      labels = c("No", "Yes")
    )
  )

str(bmt_clean)
summary(bmt_clean)
glimpse(bmt_clean)
summary(bmt_clean[c("t1", "t2", "ta", "tc", "tp")])
lapply(bmt_clean[c("d1", "d2", "d3", "da", "dc", "dp")], table)
