#### polyxsolve function ####
# model = model name
# type  = string with the two options "second" and "third", meaning second and third order polynomial
# length = integer of sample size
# y = vector of sample fluorescence data uncorrected for fluophore autofluorescence
# control = fluophore autofluorescence, so y - control = corrected sample fluorescence
# output = string with the two options "prediction" and "x", meaning response or explanatory output

polyxsolve <- function(model, type, y, length, output){ 
  
  n <- 1:length
  
  if(type == "second"){
    
    rf2 <- Vectorize(function(x){
      polyroot(c(coef(model)[1] - y[x], 
                 coef(model)[2],
                 coef(model)[3]))[1]
    })
    
    root <- as.numeric(rf2(n))
    
  } else if(type == "third"){
    
    rf3 <- Vectorize(function(x){
      polyroot(c(coef(model)[1] - y[x], 
                 coef(model)[2],
                 coef(model)[3],
                 coef(model)[4]))[1]
    })
    
    root <- as.numeric(rf3(n))
    
  } else {
    "Not a valid model type"
  }

  prediction <- as.numeric(predict(model, newdata = data.frame(S = root)))
  
  if(output == "x"){
    root
  } else if(output == "prediction"){
    prediction
  } else {
    "Not a valid output."
  }
}

#### 220630_LAMC ####
# uncorrected for quenching
standard <- read.delim(pipe("pbpaste")) # load standard curve data from clipboard
m1 <- with(standard, lm(fluo ~ poly(S, 2, raw = T)))
coef(m1)
summary(m1)

y <- read.delim(pipe("pbpaste"), header = F)$V1 # load sample data from clipboard

# test accuracy of function output
polyxsolve(m1, "second", y, 3, "prediction")
y 

# find x values
x <- polyxsolve(m1, "second", y, 3, "x")
x

# export x values
write.csv(x, "~/Desktop/x.csv")

# corrected for quenching
standard <- read.delim(pipe("pbpaste")) # load standard curve data from clipboard
m1 <- with(standard, lm(fluo ~ poly(S, 2, raw = T)))
coef(m1)
summary(m1)

y <- read.delim(pipe("pbpaste"), header = F)$V1 # load sample data from clipboard

# test accuracy of function output
polyxsolve(m1, "second", y, 45, "prediction")
y 

# find x values
x <- polyxsolve(m1, "second", y, 45, "x")
x

# export x values
write.csv(x, "~/Desktop/x.csv")


#### 220706_LAMC ####
# uncorrected for quenching
standard <- read.delim(pipe("pbpaste")) # load standard curve data from clipboard
m1 <- with(standard, lm(fluo ~ poly(S, 2, raw = T)))
coef(m1)
summary(m1)

y <- read.delim(pipe("pbpaste"), header = F)$V1 # load sample data from clipboard

# test accuracy of function output
polyxsolve(m1, "second", y, 3, "prediction")
y 

# find x values
x <- polyxsolve(m1, "second", y, 3, "x")
x

# export x values
write.csv(x, "~/Desktop/x.csv")

# corrected for quenching
standard <- read.delim(pipe("pbpaste")) # load standard curve data from clipboard
m1 <- with(standard, lm(fluo ~ poly(S, 2, raw = T)))
coef(m1)
summary(m1)

y <- read.delim(pipe("pbpaste"), header = F)$V1 # load sample data from clipboard

# test accuracy of function output
polyxsolve(m1, "second", y, 48, "prediction")
y 

# find x values
x <- polyxsolve(m1, "second", y, 48, "x")
x

# export x values
write.csv(x, "~/Desktop/x.csv")



#### 220713_LAMC ####
# uncorrected for quenching
standard <- read.delim(pipe("pbpaste")) # load standard curve data from clipboard
m1 <- with(standard, lm(fluo ~ poly(S, 2, raw = T)))
coef(m1)
summary(m1)

y <- read.delim(pipe("pbpaste"), header = F)$V1 # load sample data from clipboard

# test accuracy of function output
polyxsolve(m1, "second", y, 3, "prediction")
y 

# find x values
x <- polyxsolve(m1, "second", y, 3, "x")
x

# export x values
write.csv(x, "~/Desktop/x.csv")

# corrected for quenching
standard <- read.delim(pipe("pbpaste")) # load standard curve data from clipboard
m1 <- with(standard, lm(fluo ~ poly(S, 2, raw = T)))
coef(m1)
summary(m1)

y <- read.delim(pipe("pbpaste"), header = F)$V1 # load sample data from clipboard

# test accuracy of function output
polyxsolve(m1, "second", y, 48, "prediction")
y 

# find x values
x <- polyxsolve(m1, "second", y, 48, "x")
x

# export x values
write.csv(x, "~/Desktop/x.csv")


#### 220727_LAMC ####
# uncorrected for quenching
standard <- read.delim(pipe("pbpaste")) # load standard curve data from clipboard
m1 <- with(standard, lm(fluo ~ poly(S, 2, raw = T)))
coef(m1)
summary(m1)

y <- read.delim(pipe("pbpaste"), header = F)$V1 # load sample data from clipboard

# test accuracy of function output
polyxsolve(m1, "second", y, 3, "prediction")
y 

# find x values
x <- polyxsolve(m1, "second", y, 3, "x")
x

# export x values
write.csv(x, "~/Desktop/x.csv")

# corrected for quenching
standard <- read.delim(pipe("pbpaste")) # load standard curve data from clipboard
m1 <- with(standard, lm(fluo ~ poly(S, 2, raw = T)))
coef(m1)
summary(m1)

y <- read.delim(pipe("pbpaste"), header = F)$V1 # load sample data from clipboard

# test accuracy of function output
polyxsolve(m1, "second", y, 39, "prediction")
y 

# find x values
x <- polyxsolve(m1, "second", y, 39, "x")
x

# export x values
write.csv(x, "~/Desktop/x.csv")



#### 220806_LAMC ####
# uncorrected for quenching
standard <- read.delim(pipe("pbpaste")) # load standard curve data from clipboard
m1 <- with(standard, lm(fluo ~ poly(S, 2, raw = T)))
coef(m1)
summary(m1)

y <- read.delim(pipe("pbpaste"), header = F)$V1 # load sample data from clipboard

# test accuracy of function output
polyxsolve(m1, "second", y, 3, "prediction")
y 

# find x values
x <- polyxsolve(m1, "second", y, 3, "x")
x

# export x values
write.csv(x, "~/Desktop/x.csv")

# corrected for quenching
standard <- read.delim(pipe("pbpaste")) # load standard curve data from clipboard
m1 <- with(standard, lm(fluo ~ poly(S, 2, raw = T)))
coef(m1)
summary(m1)

y <- read.delim(pipe("pbpaste"), header = F)$V1 # load sample data from clipboard

# test accuracy of function output
polyxsolve(m1, "second", y, 27, "prediction")
y 

# find x values
x <- polyxsolve(m1, "second", y, 27, "x")
x

# export x values
write.csv(x, "~/Desktop/x.csv")


#### 221026_LAMC ####
# uncorrected for quenching
standard <- read.delim(pipe("pbpaste")) # load standard curve data from clipboard
m1 <- with(standard, lm(fluo ~ S))
coef(m1)
summary(m1)

# corrected for quenching
standard <- read.delim(pipe("pbpaste")) # load standard curve data from clipboard
m1 <- with(standard, lm(fluo ~ poly(S, 2, raw = T)))
coef(m1)
summary(m1)

y <- read.delim(pipe("pbpaste"), header = F)$V1 # load sample data from clipboard

# test accuracy of function output
polyxsolve(m1, "second", y, 6, "prediction")
y 

# find x values
x <- polyxsolve(m1, "second", y, 6, "x")
x

# export x values
write.csv(x, "~/Desktop/x.csv")



#### 220630_4MUG ####
# uncorrected for quenching
standard <- read.delim(pipe("pbpaste")) # load standard curve data from clipboard
m2 <- with(standard, lm(fluo ~ poly(S, 2, raw = T)))
coef(m2)
summary(m2)

y <- read.delim(pipe("pbpaste"), header = F)$V1 # load sample data from clipboard

# test accuracy of function output
polyxsolve(m2, "second", y, 3, "prediction")
y 

# find x values
x <- polyxsolve(m2, "second", y, 3, "x")
x

# export x values
write.csv(x, "~/Desktop/x.csv")


# corrected for quenching
standard <- read.delim(pipe("pbpaste")) # load standard curve data from clipboard
m3 <- with(standard, lm(fluo ~ poly(S, 2, raw = T)))
coef(m3)
summary(m3)

y <- read.delim(pipe("pbpaste"), header = F)$V1 # load sample data from clipboard

# test accuracy of function output
polyxsolve(m3, "second", y, 45, "prediction")
y 

# find x values
x <- polyxsolve(m3, "second", y, 45, "x")
x

# export x values
write.csv(x, "~/Desktop/x.csv")



#### 220713_4MUG ####
# uncorrected for quenching
standard <- read.delim(pipe("pbpaste")) # load standard curve data from clipboard
m2 <- with(standard, lm(fluo ~ poly(S, 2, raw = T)))
coef(m2)
summary(m2)

y <- read.delim(pipe("pbpaste"), header = F)$V1 # load sample data from clipboard

# test accuracy of function output
polyxsolve(m2, "second", y, 3, "prediction")
y 

# find x values
x <- polyxsolve(m2, "second", y, 3, "x")
x

# export x values
write.csv(x, "~/Desktop/x.csv")


# corrected for quenching
standard <- read.delim(pipe("pbpaste")) # load standard curve data from clipboard
m3 <- with(standard, lm(fluo ~ poly(S, 2, raw = T)))
coef(m3)
summary(m3)

y <- read.delim(pipe("pbpaste"), header = F)$V1 # load sample data from clipboard

# test accuracy of function output
polyxsolve(m3, "second", y, 48, "prediction")
y 

# find x values
x <- polyxsolve(m3, "second", y, 48, "x")
x

# export x values
write.csv(x, "~/Desktop/x.csv")


#### 220727_4MUG ####
# uncorrected for quenching
standard <- read.delim(pipe("pbpaste")) # load standard curve data from clipboard
m4 <- with(standard, lm(fluo ~ poly(S, 2, raw = T)))
coef(m4)
summary(m4)

y <- read.delim(pipe("pbpaste"), header = F)$V1 # load sample data from clipboard

# test accuracy of function output
polyxsolve(m4, "second", y, 3, "prediction")
y 

# find x values
x <- polyxsolve(m4, "second", y, 3, "x")
x

# export x values
write.csv(x, "~/Desktop/x.csv")


# corrected for quenching
standard <- read.delim(pipe("pbpaste")) # load standard curve data from clipboard
m5 <- with(standard, lm(fluo ~ poly(S, 2, raw = T)))
coef(m5)
summary(m5)

y <- read.delim(pipe("pbpaste"), header = F)$V1 # load sample data from clipboard

# test accuracy of function output
polyxsolve(m5, "second", y, 39, "prediction")
y 

# find x values
x <- polyxsolve(m5, "second", y, 39, "x")
x

# export x values
write.csv(x, "~/Desktop/x.csv")


#### 220806_4MUG ####
# uncorrected for quenching
standard <- read.delim(pipe("pbpaste")) # load standard curve data from clipboard
m6 <- with(standard, lm(fluo ~ poly(S, 2, raw = T)))
coef(m6)
summary(m6)

y <- read.delim(pipe("pbpaste"), header = F)$V1 # load sample data from clipboard

# test accuracy of function output
polyxsolve(m6, "second", y, 3, "prediction")
y 

# find x values
x <- polyxsolve(m6, "second", y, 3, "x")
x

# export x values
write.csv(x, "~/Desktop/x.csv")


# corrected for quenching
standard <- read.delim(pipe("pbpaste")) # load standard curve data from clipboard
m7 <- with(standard, lm(fluo ~ poly(S, 2, raw = T)))
coef(m7)
summary(m7)

y <- read.delim(pipe("pbpaste"), header = F)$V1 # load sample data from clipboard

# test accuracy of function output
polyxsolve(m7, "second", y, 27, "prediction")
y 

# find x values
x <- polyxsolve(m7, "second", y, 27, "x")
x

# export x values
write.csv(x, "~/Desktop/x.csv")


#### 221026_4MUG ####
# uncorrected for quenching
standard <- read.delim(pipe("pbpaste")) # load standard curve data from clipboard
m6 <- with(standard, lm(fluo ~ poly(S, 2, raw = T)))
coef(m6)
summary(m6)

y <- read.delim(pipe("pbpaste"), header = F)$V1 # load sample data from clipboard

# test accuracy of function output
polyxsolve(m6, "second", y, 3, "prediction")
y 

# find x values
x <- polyxsolve(m6, "second", y, 3, "x")
x

# export x values
write.csv(x, "~/Desktop/x.csv")


# corrected for quenching
standard <- read.delim(pipe("pbpaste")) # load standard curve data from clipboard
m7 <- with(standard, lm(fluo ~ poly(S, 2, raw = T)))
coef(m7)
summary(m7)

y <- read.delim(pipe("pbpaste"), header = F)$V1 # load sample data from clipboard

# test accuracy of function output
polyxsolve(m7, "second", y, 6, "prediction")
y 

# find x values
x <- polyxsolve(m7, "second", y, 6, "x")
x

# export x values
write.csv(x, "~/Desktop/x.csv")
