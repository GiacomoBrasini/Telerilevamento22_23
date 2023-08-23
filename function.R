# FUNCTIONS in R

# greet
greet <- function() {
  print("Hello")
}
greet()

# greet me
greet_me <- function(name) {
  string <- paste("Hello", name, sep = " ")
  return(string)
}
greet_me("Giacomo")

# greet me with default
greet_me <- function(name = "Duccio") {
  string <- paste("Hello", name, sep = " ")
  return(string)
}
greet_me() # Duccio is default, but it can be changed
greet_me("Giacomo")

# add two to a number
add_two <- function(num = 0) {
  n_plus2 <- num + 2
  return(n_plus2)
}
add_two(3)
add_two("f")

# if else 
n <- 3
if(n < 10) {
  print("n is less than 10")
} else{
  print("n is greater than or equal to 10")
}

n <- 10
if(n < 10) {
  print("n is less than 10")
} else if(n >10) {
  print("n is greater than 10")
} else{
  print("n is equal to 10")
}

n <- 33
if(n < 10) {
  print("n is less than 10")
} else if(n >10) {
  print("n is greater than 10")
} else{
  print("n is equal to 10")
}

# We can create a function that tell if a choosen number is >, < or = to 10
is_ten <- function(num = 0){
  if(num < 10) {
    print("n is less than 10")
  } else if(num >100) {
    print("n is greater than 10")
  } else{
    print("n is equal to 10")
  }
}
is_ten(121)
is_ten(4)
is_ten(10)

# add 2 function
add_two <- function(n = 0) {
  if(is.numeric(n)){
    n_plus2 <- n + 2
    return(n_plus2)
  } else if(!is.numeric(n)){
    print("n is not a number")
  }
}
add_two(3)
add_two("f")
add_two(373)

# another way to do it, using stop to report an error
add_two <- function(n = 0) {
  if(!is.numeric(n)) stop("n is not numeric")
  n_plus2 <- n + 2
  return(n_plus2)
}
add_two(3)
add_two("f")

# we can report an NA using warning
add_two <- function(n = 0) {
  if(is.na(n)) warning("n is NA")
  if(!is.numeric(n)) stop("n is not numeric")
  n_plus2 <- n + 2
  return(n_plus2)
}
add_two(3)
add_two("f")
add_two(NA)

# functions create an isolated environment


# For loops

for (i in 1:5) {
  print(i)
}

num <- 0
for (i in 1:3) {
  num <- num + 5
  print(num)
}

num <- 1
for (i in 1:4) {
  num <- num + 5
  print(num)
}

# Function including a for loop
multiply <- function(n1 = 1, n2 = 1){
  if(!is.numeric(n1) || !is.numeric(n2)) stop("n1 or n2 or both are not numeric")
  num <- 0
  for (i in 1:n1) {
    num <- num + n2
  }
  
  return(num)
}

multiply(5, 12)
