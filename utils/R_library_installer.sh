#!/usr/bin/env Rscript
# pluta 10/25/12

# run this script first to install all necessary R libraries.



pckg = try(require(Rcpp))

if(!pckg)
{
	cat("Installing 'Rcpp' from CRAN\n")
	install.packages("Rcpp", repos = "http://cran.r-project.org")
	require("Rcpp")
}
	
	
pckg = try(require(signal))
if(!pckg)
{
	cat("Installing 'signal' from CRAN\n")
	install.packages("signal", repos = "http://cran.r-project.org")
	require("signal")
}		

pckg = try(require(timeSeries))
if(!pckg)
{
	cat("Installing 'timeSeries' from CRAN\n")
	install.packages("timeSeries", repos = "http://cran.r-project.org")
	require("timeSeries")
}	

pckg = try(require(mFilter))
if(!pckg)
{
	cat("Installing 'mFilter' from CRAN\n")
	install.packages("mFilter", repos = "http://cran.r-project.org")
	require("mFilter")
}	

pckg = try(require(doParallel))
if(!pckg)
{
	cat("Installing 'doParallel' from CRAN\n")
	install.packages("doParallel", repos = "http://cran.r-project.org")
	require("doParallel")
}	

pckg = try(require(robust))
if(!pckg)
{
	cat("Installing 'robust' from CRAN\n")
	install.packages("robust", repos = "http://cran.r-project.org")
	require("robust")
}	

pckg = try(require(brainwaver))
if(!pckg)
{
	cat("Installing 'brainwaver' from CRAN\n")
	install.packages("brainwaver", repos = "http://cran.r-project.org")
	require("brainwaver")
}	

pckg = try(require(igraph))
if(!pckg)
{
	cat("Installing 'igraph' from CRAN\n")
	install.packages("igraph", repos = "http://cran.r-project.org")
	require("igraph")
}	