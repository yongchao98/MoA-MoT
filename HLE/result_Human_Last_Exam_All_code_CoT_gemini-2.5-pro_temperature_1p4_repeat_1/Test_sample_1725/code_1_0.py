import numpy as np

# Part A: Dimension of the fibers of the bundle pi
# The connection C is a general affine connection on R^3. Its components are C^i_{jk}.
# To describe dislocations via torsion, the connection cannot be assumed symmetric in the lower two indices.
# Thus, i, j, k each run from 1 to 3, giving 3*3*3 components.
dim_fiber = 3 * 3 * 3

# Part B: Number of coefficients specifying the energy E
# E is a quadratic form on the 9-dimensional space of the torsion tensor T.
# The number of coefficients is the number of invariants in the symmetric square of the
# representation of the torsion tensor under the cubic group O_h.
# This can be calculated using character theory.
# Number of invariants = (1/|G|) * sum_over_classes(class_order * character)
# where the character is for the symmetric square of the representation.

# Group O_h has order 48.
group_order = 48

# Class orders for O_h: E, 8C_3, 6C_2, 6C_4, 3C_2', i, 8S_6, 6sigma_d, 6S_4, 3sigma_h
# Note: There are variations in naming/ordering of classes in different character tables.
# The calculation below uses a consistent table. The final result is robust to these choices.
class_orders = np.array([1, 8, 6, 6, 3, 1, 8, 6, 6, 3])

# The character of the symmetric square of the torsion representation (T_1u x T_1g)
# has been calculated to be the following for the 10 classes of O_h.
# Chi_Sym2(g) = 0.5 * ( Chi(g)^2 + Chi(g^2) )
# The result of this calculation is:
chi_sym2 = np.array([45, 0, 5, 1, 5, 45, 0, 5, 1, 5])

# The number of invariants is the inner product of the character vector with the
# character of the trivial representation A_1g (which is all 1s), divided by the group order.
num_invariants = (1 / group_order) * np.sum(class_orders * chi_sym2)

# We print the final answer in the required format.
# The final values for the equation are printed to show the calculation.
print("Calculation for Part B:")
calculation_str = f"(1/{group_order}) * ("
for i in range(len(class_orders)):
    calculation_str += f"{class_orders[i]}*{chi_sym2[i]}"
    if i < len(class_orders) - 1:
        calculation_str += " + "
calculation_str += f") = {int(num_invariants)}"
print(calculation_str)

print("\nFinal Answer:")
# The format should be "A B"
print(f"{dim_fiber} {int(num_invariants)}")
<<<27 4>>>