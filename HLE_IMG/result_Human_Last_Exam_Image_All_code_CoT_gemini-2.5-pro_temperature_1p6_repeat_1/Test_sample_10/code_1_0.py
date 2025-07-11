#
# This script predicts the inversion barrier for a PAH molecule based on a linear trend.
#

# --- Step 1: Define the known data points ---
# Let n be the number of 5-membered rings in the molecule.
# Let B be the inversion barrier in kcal/mol.

# Data for the first molecule (dibenzo[ghi,mno]fluoranthene)
n1 = 1
B1 = 10

# Data for the second molecule (diacenaphtho[...]chrysene)
n2 = 2
B2 = 49

# The number of 5-membered rings for the target molecule (triacenaphtho[...]triphenylene)
n3 = 3

# --- Step 2: Assume a linear relationship B(n) = a*n + b and find a and b ---
# We have a system of two linear equations:
# a*1 + b = 10
# a*2 + b = 49

# The slope 'a' is the change in B divided by the change in n.
a = (B2 - B1) / (n2 - n1)

# The intercept 'b' can be found using the first data point.
# b = B1 - a*n1
b = B1 - a * n1

# --- Step 3: Calculate the predicted barrier for the target molecule (n=3) ---
predicted_barrier = a * n3 + b

# --- Step 4: Print the result including the final equation ---
# The calculated coefficients 'a' and 'b' are integers.
# The predicted barrier is also an integer, as requested.
print("Based on a linear extrapolation from the provided data:")
print(f"The predicted inversion barrier is calculated using the equation: {int(a)} * n - {abs(int(b))}")
print("For the target molecule with n = 3:")
print(f"Predicted Barrier = {int(a)} * {n3} - {abs(int(b))} = {int(predicted_barrier)} kcal/mol")