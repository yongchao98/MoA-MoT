import math

# Given values
T = math.log(10)
# We calculate B from the expression, although we will use the simplified analytic form later
# B = (1/2) * (10**20) / (99**2)

# From the derivation, we have the equation for alpha^8:
# alpha^8 = (99^6 * 10^4) / 64

# Calculate alpha^8
alpha_8 = (99**6 * 10**4) / 64

# Calculate alpha by taking the eighth root
alpha = alpha_8**(1/8)

# The result is alpha = (99^(6/8) * 10^(4/8)) / (64^(1/8))
# alpha = (99^(3/4) * 10^(1/2)) / ( (2^6)^(1/8) )
# alpha = (99**(3/4) * sqrt(10)) / (2**(6/8))
# alpha = (99**(3/4) * sqrt(10)) / (2**(3/4))
# alpha = (99/2)**(3/4) * sqrt(10)

print(f"The equation for alpha^8 is: alpha^8 = (99^6 * 10^4) / 64")
print(f"99^6 = {99**6}")
print(f"10^4 = {10**4}")
print(f"99^6 * 10^4 = {99**6 * 10**4}")
print(f"alpha^8 = {99**6 * 10**4} / 64 = {alpha_8}")
print(f"alpha = ({alpha_8})^(1/8) = {alpha}")

# Final Answer
# print("<<<{}>>>".format(alpha))
<<<59.01344243169823>>>