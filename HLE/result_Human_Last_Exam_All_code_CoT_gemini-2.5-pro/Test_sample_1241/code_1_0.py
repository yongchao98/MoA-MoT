import sys

# Step 1: Define the lambda values
lambda01 = 0.019
lambda10 = 0.65
lambda12 = 0.4
lambda21 = 0.392
lambda23 = 0.008
lambda31 = 0.008

# Step 2: Solve the system of linear equations for the steady-state probabilities
# We express P0, P2, P3 in terms of P1.
# Let P_i = c_i * P1
# P1 = c1 * P1 => c1 = 1
c1 = 1.0

# From l01*P0 = l10*P1
c0 = lambda10 / lambda01

# From (l21 + l23)*P2 = l12*P1
c2 = lambda12 / (lambda21 + lambda23)

# From l31*P3 = l23*P2
c3 = (lambda23 / lambda31) * c2

# The normalization condition is P0 + P1 + P2 + P3 = 1
# (c0*P1) + (c1*P1) + (c2*P1) + (c3*P1) = 1
# P1 * (c0 + c1 + c2 + c3) = 1
# P1 = 1 / (c0 + c1 + c2 + c3)
sum_of_coeffs = c0 + c1 + c2 + c3
P1 = 1 / sum_of_coeffs

# Now calculate P0
P0 = c0 * P1

# Step 3: Calculate the required sum P0 + P1
result = P0 + P1

# Print the final equation with the calculated values
# Using f-strings to format the output with high precision
if sys.version_info >= (3, 6):
    print(f"P0(+inf) + P1(+inf) = {P0:.9f} + {P1:.9f} = {result:.9f}")
else:
    # Fallback for older Python versions
    print("P0(+inf) + P1(+inf) = {:.9f} + {:.9f} = {:.9f}".format(P0, P1, result))
