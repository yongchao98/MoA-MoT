import math

# Step 1 & 2: Define the relationship between p, n, and the minimal decay exponent alpha_min.
# The function I(a) is not in L^p(R^n) for p <= n / alpha_min.
# We want the largest such p, which is p = n / alpha_min.

# Step 3: Determine the dimension of the parameter space, n.
# The integral I depends on 9 parameters, a_1 through a_9.
n = 9

# Step 4, 5, & 6: Determine the minimal decay exponent, alpha_min.
# The decay of the integral I(a) is slowest when the phase polynomial is as flat as possible.
# This occurs when the phase reduces to a single variable with the highest possible degeneracy.
# For a phase like a_9 * y^3, the integral decays as |a_9|^(-1/3).
# The highest order of the polynomial is 3, so the minimal decay exponent is 1/3.
alpha_min_num = 1
alpha_min_den = 3
alpha_min = alpha_min_num / alpha_min_den

# Step 7: Calculate the final value of p.
p = n / alpha_min

print("The final equation for p is:")
print(f"p = n / alpha_min = {n} / ({alpha_min_num}/{alpha_min_den}) = {int(p)}")