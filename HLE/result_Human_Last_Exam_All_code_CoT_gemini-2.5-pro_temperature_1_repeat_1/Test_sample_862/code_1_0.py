import numpy as np

# Step 1: Define the equation to solve for k = sqrt(lambda).
# The equation is tan(k * pi / 4) = 1 / sqrt(3).
# We print the components of this equation.
val_tan = 1 / np.sqrt(3)
a = 4
b = 1
c = 3
print(f"The transcendental equation to solve for k=sqrt(lambda) is tan(k*pi/{a}) = {b}/sqrt({c})")
print(f"The value of the right hand side is {val_tan}")

# Step 2: Solve for the angle alpha = k * pi / 4.
# The inverse tangent gives the principal value.
alpha = np.arctan(val_tan)
print(f"The angle alpha = k*pi/4 in radians is {alpha:.4f}, which is pi/6.")

# Step 3: Solve for k.
# k = alpha * 4 / pi
k = (np.pi/6) * 4 / np.pi
numerator_k = 2
denominator_k = 3
print(f"Solving for k gives k = {numerator_k}/{denominator_k}")
print(f"The value of k = sqrt(lambda) is {k}")

# Step 4: Calculate lambda.
# lambda = k^2
lambda_1 = k**2
numerator_lambda = 4
denominator_lambda = 9
print(f"The eigenvalue lambda is k^2 = ({numerator_k}/{denominator_k})^2 = {numerator_lambda}/{denominator_lambda}")
print(f"The value of lambda is {lambda_1}")

# Step 5: Calculate the constant C.
# C = 1 / lambda_1
C = 1 / lambda_1
numerator_C = 9
denominator_C = 4
print(f"The constant C is 1/lambda = 1/({numerator_lambda}/{denominator_lambda}) = {numerator_C}/{denominator_C}")
print(f"The final value of C is {C}")
