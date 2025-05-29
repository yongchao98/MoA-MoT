from sympy import sympify, degree_list, symbols

# Define the equation
equation_str = "x**2 + y + z"

# Convert the input string to a SymPy expression
n1 = sympify(equation_str)

# Get the degrees of the variables in the equation
x, y, z = symbols('x y z')
degrees = [n1.degree(x), n1.degree(y), n1.degree(z)]

# Determine the type of equation based on the degrees
if max(degrees) == 2:
    result = "quadratic equation"
elif max(degrees) == 1:
    if degrees.count(1) == 3:
        result = "linear equation of three variables"
    elif degrees.count(1) == 2:
        result = "linear equation of two variables"
    else:
        result = "linear equation"
elif max(degrees) == 3:
    result = "cubic equation"
elif max(degrees) == 4:
    result = "biquadratic equation"
elif max(degrees) > 4:
    result = "polynomial equation of degree greater than 4"
else:
    result = "unknown equation type"

print(result)