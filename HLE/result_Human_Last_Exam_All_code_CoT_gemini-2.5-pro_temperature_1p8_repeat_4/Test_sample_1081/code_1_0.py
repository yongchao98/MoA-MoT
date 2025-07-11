# This script calculates the number of F_q-rational maximal tori classes
# for a reductive group G of type E_8 over the finite field F_q.

# The number of such classes is equal to the number of conjugacy classes
# of the Weyl group of type E_8, denoted W(E_8).
num_conjugacy_classes_W_E8 = 112

# The problem is interpreted as finding the number of G(F_q)-conjugacy classes of
# F_q-rational maximal tori. This number is constant and independent of q.

# The final equation is derived from the fact that for a split group like E_8,
# the number of tori classes equals the number of Weyl group conjugacy classes.
#
# Equation:
# Number of rational maximal tori classes = Number of W(E_8) conjugacy classes
# Number of rational maximal tori classes = 112

print("The final equation for the number of classes of rational maximal tori is:")

# We output the final equation with the computed number.
print(f"Number of rational maximal tori classes = {num_conjugacy_classes_W_E8}")
