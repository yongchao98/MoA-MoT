import math

# This script calculates the leading-order asymptotic mass ratio between the
# lightest and subsequent solitonic excitations in the CP(N-1) model.

# Step 1: The mass M_k of the k-th solitonic excitation is proportional to k/sin(pi*k/N).
# We are interested in the ratio for the lightest (k=1) and subsequent (k=2) excitations.
k_lightest = 1
k_subsequent = 2

# Step 2: The mass ratio R(N) is given by:
# R(N) = (k_subsequent / sin(pi * k_subsequent / N)) / (k_lightest / sin(pi * k_lightest / N))
# After substituting k=1 and k=2 and simplifying, we use the identity sin(2x) = 2*sin(x)*cos(x).
# This gives: R(N) = 1 / cos(pi / N)

# Step 3: We need to find the leading-order asymptotic ratio as N approaches infinity.
# This means we take the limit of R(N) as N -> infinity.
# In this limit, the argument of cosine, (pi / N), approaches 0.

# Step 4: Calculate the values needed for the final equation.
# The numerator in the simplified ratio is a constant.
numerator = 1

# The denominator in the limit is cos(0).
# math.cos(0) gives the value of cos(0).
limit_denominator = math.cos(0)

# The asymptotic ratio is the numerator divided by the limit of the denominator.
asymptotic_ratio = numerator / limit_denominator

# Step 5: Print the components of the final calculation and the result.
# The final equation is: Ratio = 1 / 1 = 1.
print("The leading-order asymptotic mass ratio is calculated from the limit of 1/cos(pi/N) as N -> infinity.")
print("The final equation for the ratio is derived as follows:")
print(f"Numerator = {numerator}")
print(f"Denominator = cos(0) = {int(limit_denominator)}")
print("Resulting Equation:")
print(f"{numerator} / {int(limit_denominator)} = {int(asymptotic_ratio)}")
