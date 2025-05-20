#' ---
#' title: "Intro to R"
#' author: "Hyejeong Choi"
#' date: "April 10, 2025"
#' ---
# My first R script
x <- 1:50
plot(x, sin(x))
plot(x, sin(x), typ="l", col="green", lwd=3, 
     xlab="Silly x axis", ylab="Sensible y axis")