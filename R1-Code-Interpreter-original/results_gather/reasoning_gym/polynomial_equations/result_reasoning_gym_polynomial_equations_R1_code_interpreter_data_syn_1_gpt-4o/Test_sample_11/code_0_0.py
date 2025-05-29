import math

# Coefficients
A = 86
B = -2
C = -13

# Calculate the discriminant
D = B**2 - 4*A*C

# Check the discriminant and calculate the solutions
if D > 0:
    # Two real solutions
    b1 = (-B + math.sqrt(D)) / (2*A)
    b2 = (-B - math.sqrt(D)) / (2*A)
    result = f"{b1:.4f}, {b2:.4f}"
elif D == 0:
    # One real solution
    b = -B / (2*A)
    result = f"{b:.4f}"
else:
    # No real solutions
    result = ""

print(result)