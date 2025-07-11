# Based on the set-theoretic analysis, we determined the following:
# delta_1 = supremum(X) = omega_2
# delta_2 = infimum(X) = omega_2

# We can represent the cardinal omega_n by its index n.
# So, omega_2 is represented by the integer 2.
delta_1_index = 2
delta_2_index = 2

# Cardinal addition for infinite cardinals is defined as the maximum.
# So, omega_a + omega_b = omega_{max(a, b)}.
result_index = max(delta_1_index, delta_2_index)

# The problem asks for the sum delta_1 + delta_2. The result is omega_2.
# The final equation is omega_2 + omega_2 = omega_2.
# We print each cardinal that appears in this final equation.
print(f"omega_{delta_1_index} + omega_{delta_2_index} = omega_{result_index}")