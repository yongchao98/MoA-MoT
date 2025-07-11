import math

# Based on the asymptotic analysis of the ODE for large x, the solution
# takes the form y(x) = A * x^n.
# The analysis yields n = -1 and an equation for A, which is A^3 = -6.

# 1. Define the equation for the coefficient A.
A_cubed = -6.0

# 2. Solve for A by calculating the real cube root of -6.
# The sign of A will be negative.
A = math.copysign(pow(abs(A_cubed), 1.0/3.0), A_cubed)

# 3. Round the value of A to two decimal places as requested.
A_rounded = round(A, 2)

# 4. The power n was determined to be -1 from the dominant balance.
n = -1

# 5. Construct and print the final analytical expression.
# The format y(x) = C * x^(n) explicitly shows both numbers in the equation.
print(f"y(x) = {A_rounded} * x^({n})")