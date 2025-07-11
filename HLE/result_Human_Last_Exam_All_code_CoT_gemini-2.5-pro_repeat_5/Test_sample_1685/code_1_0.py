import math

# Step 1: Define the constant c from the theoretical derivation c^3 = 10/9
c = (10/9)**(1/3)

# Step 2: Calculate the singularity location x_0 using the initial condition y'(0)=3
# x_0 = (c/3)^(3/2)
x_0 = (10/243)**(1/2) # Simplified exact form

# Step 3: Calculate the integration constant K using the initial condition y(0)=0
# K = 3*c*x_0**(1/3)
K = math.sqrt(10/3) # Simplified exact form

# Step 4: The final expression is y(x) = K - 3*c*(x_0 - x)^(1/3)
# We will print the equation with the coefficients rounded to two decimal places.

coeff_K = K
coeff_3c = 3 * c
coeff_x0 = x_0

print("The analytical expression that approximates the solution is:")
# The final equation has the form y(x) = A - B * (C - x)^(1/3)
# We output each number in this final equation format.
print(f"y(x) = {coeff_K:.2f} - {coeff_3c:.2f} * ({coeff_x0:.2f} - x)^(1/3)")
