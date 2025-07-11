# Step 1: Determine the constants K1, K2, K3, K4
# The denominator of the coefficient A_m is calculated by the integral of sin^2(sqrt(lambda_m)*x) from 0 to l.
# Calculated denominator = l/2 - sin(2*l*sqrt(lambda_m))/(4*sqrt(lambda_m))
# Given denominator form = (1/(K1*sqrt(lambda_m))) * (K2*l*sqrt(lambda_m) + K3*sin(K4*l*sqrt(lambda_m)))
#
# By comparing the calculated denominator with the given formula, we equate them:
# l/2 - sin(2*l*sqrt(lambda_m))/(4*sqrt(lambda_m)) = (K2*l)/K1 + (K3*sin(K4*l*sqrt(lambda_m)))/(K1*sqrt(lambda_m))
#
# From the argument of the sine function, we deduce K4:
# K4 * l * sqrt(lambda_m) = 2 * l * sqrt(lambda_m) => K4 = 2
K4 = 2
#
# By comparing the other terms, we find relationships between K1, K2, and K3:
# For the term independent of sin(): l/2 = (K2*l)/K1 => K1 = 2*K2
# For the sin() term: -1/(4*sqrt(lambda_m)) = K3/(K1*sqrt(lambda_m)) => K1 = -4*K3
#
# These relations do not uniquely determine the constants, but their ratios are fixed.
# A standard practice in such problems is to choose the simplest set of integers.
# Let's set the scaling by choosing K3 = -1.
K3 = -1
#
# Now we can find K1 and K2.
# From K1 = -4*K3:
K1 = -4 * K3
# From K1 = 2*K2:
K2 = K1 / 2
# Let's cast K2 to an integer as it should be.
K2 = int(K2)

# Step 2: Determine the constant K.
# The value of K is not explicitly defined in the problem, which suggests a possible typo.
# A reasonable assumption is that K is related to the ratios of the other constants.
# Let's assume K is the ratio of K3 to K2.
K = K3 / K2

# Step 3: Calculate the final product K * K1 * K2 * K3 * K4
product = K * K1 * K2 * K3 * K4

# Step 4: Print the final result in an equation format
print(f"Based on the derivation, the determined constants are:")
print(f"K1 = {K1}")
print(f"K2 = {K2}")
print(f"K3 = {K3}")
print(f"K4 = {K4}")
print(f"The value for K is inferred as K = K3/K2 = {K3}/{K2} = {K}")
print("\nThe final product is calculated as:")
print(f"Product = K * K1 * K2 * K3 * K4")
print(f"Product = ({K}) * ({K1}) * ({K2}) * ({K3}) * ({K4}) = {product}")