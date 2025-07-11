import sys

# Step 1: Define the given values.
# Firing strength of the rule's antecedent part.
phi_k = 0.7
# Membership grade of a point in the original consequent set.
mu_G = 0.9

# Step 2: Explain the method.
# A t-norm is used to combine these values. The minimum t-norm is a common choice.
# This calculates the membership grade in the final inferred fuzzy set (clipping method).
print(f"Calculating the rule activation level using the minimum t-norm.")
print(f"The inputs are the antecedent firing strength Phi_k(x′) = {phi_k} and the consequent membership grade μG_k_j(y_j) = {mu_G}.")

# Step 3: Apply the minimum t-norm.
activation_level = min(phi_k, mu_G)

# Step 4: Print the final equation and the result.
print(f"\nThe final equation is: activation_level = min({phi_k}, {mu_G})")
print(f"The resulting rule activation level is: {activation_level}")

# The following line is for the platform to extract the final answer.
# It should not be part of the user's script.
sys.stdout.write(f"\n<<<{activation_level}>>>")