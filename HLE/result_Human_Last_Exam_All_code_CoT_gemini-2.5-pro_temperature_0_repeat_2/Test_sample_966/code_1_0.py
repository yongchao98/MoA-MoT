import math

# This script calculates the dimension of the middle cohomology group H^100(X, Q)
# for a complete intersection X of two quadrics in CP^102.

# The dimension is given by b_100(X) = chi(X) - 100, where chi(X) is the Euler characteristic.
# The formula for chi(X) is: 4*C(200, 98) - 4*C(199, 97) + C(198, 96).

# Step 1: Calculate the binomial coefficients
c1 = math.comb(200, 98)
c2 = math.comb(199, 97)
c3 = math.comb(198, 96)

# Step 2: Calculate the intermediate terms for chi(X)
term1 = 4 * c1
term2 = 4 * c2
term3 = c3

# Step 3: Calculate chi(X)
chi_X = term1 - term2 + term3

# Step 4: Calculate the final dimension b_100(X)
b_100 = chi_X - 100

# Step 5: Print the final equation showing all the calculated numbers
print(f"The dimension of the middle cohomology group is given by the equation:")
print(f"dim H^100(X, Q) = (4 * C(200, 98) - 4 * C(199, 97) + C(198, 96)) - 100")
print(f"               = (4 * {c1} - 4 * {c2} + {c3}) - 100")
print(f"               = ({term1} - {term2} + {term3}) - 100")
print(f"               = {chi_X} - 100")
print(f"               = {b_100}")