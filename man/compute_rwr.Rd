% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compute_rwr.R
\name{compute_rwr}
\alias{compute_rwr}
\title{Produce rarity-weighted richness values for each site.}
\usage{
compute_rwr(data, nselect)
}
\arguments{
\item{data}{a presence/absence matrix with sites across the rows and species down the columns. Should include site names for the rows.}

\item{nselect}{number of sites to select to maximize the rarity-weighted richness score.}
}
\value{
Returns a list including a vector of the selected sites, an integer number of species included, and a vector of rarity-weighted richness scores for the selected sites.
}
\description{
Produce rarity-weighted richness values for each site.
}
\details{
This function selects sites based on rarity-weighted richness.
}
