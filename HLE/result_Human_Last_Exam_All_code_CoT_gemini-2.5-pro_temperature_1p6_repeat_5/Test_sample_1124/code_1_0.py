import math

# The value of R for which the standard map is chaotic
R = 3.57

# The modified logistic map is X_n+1 = R * X_n * (1 - X_n) + R/3
# To find the equilibrium point X, we solve the equation:
# X = R * X * (1 - X) + R/3
# This can be rearranged into the quadratic form a*X^2 + b*X + c = 0

# Coefficients of the quadratic equation
a = R
b = 1 - R
c = -R / 3

print("The modified logistic map is: X_n+1 = R * X_n * (1 - X_n) + R/3")
print(f"For R = {R}, we solve the equilibrium equation: X = {R}*X*(1-X) + {R}/3")
print("This simplifies to the following quadratic equation:")
# The problem asks to output each number in the final equation.
# The equation is aX^2 + bX + c = 0
# We will print it with the signs handled correctly.
print(f"{a:.2f}*X^2 + ({b:.2f})*X + ({c:.2f}) = 0")
print("")

# Calculate the discriminant
discriminant = b**2 - 4 * a * c

# Solve for the two possible values of X (the roots)
if discriminant >= 0:
    x1 = (-b + math.sqrt(discriminant)) / (2 * a)
    x2 = (-b - math.sqrt(discriminant)) / (2 * a)
    
    # We are looking for the equilibrium point approximately equal to 1.05
    # The positive root is the one we are interested in.
    equilibrium_point = x1
    
    print("The equilibrium point is the positive root of this equation.")
    print(f"Equilibrium Point X = {equilibrium_point}")
else:
    print("The equation has no real roots for the equilibrium point.")
