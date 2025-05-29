from sympy import sympify, degree_list, symbols

# Define the equation
equation_str = "x**2 + y + z"

# Convert the input string to a SymPy expression
n1 = sympify(equation_str)

# Get the degrees of the variables in the equation
degrees = degree_list(n1)

# Determine the type of equation based on the degrees
if len(degrees) == 1:
    p = degrees[0]
    q = 0
    r = 0
elif len(degrees) == 2:
    p = degrees[0]
    q = degrees[1]
    r = 0
elif len(degrees) == 3:
    p = degrees[0]
    q = degrees[1]
    r = degrees[2]
else:
    result = "Invalid equation"

# Determine the type of equation
if p == 1 and q == 1 and r == 1:
    result = "linear equation of three variables"
elif (p == 1 and q == 1) or (q == 1 and r == 1) or (p == 1 and r == 1):
    result = "linear equation of two variables"
else:
    if p == 1 or q == 1 or r == 1:
        result = "linear equation"
    if p == 2 or q == 2 or r == 2:
        result = "quadratic equation"
    if p == 3 or q == 3 or r == 3:
        result = "cubic equation"
    if p == 4 or q == 4 or r == 4:
        result = "biquadratic equation"
    if p > 4 or q > 4 or r > 4:
        result = "polynomial equation of degree greater than 4"

print(result)