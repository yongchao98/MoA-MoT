# This script calculates the sum of all entries in the Cartan matrix
# for the principal block of kG, where G = A_5 x C_2 and char(k) = 2.

# Based on representation theory, this sum is equal to the order of the
# Sylow 2-subgroup of G.

# The order of A_5 is 60 = 2^2 * 15.
# The order of the Sylow 2-subgroup of A_5 is 4.
order_sylow_A5 = 4

# The order of C_2 is 2.
# The order of the Sylow 2-subgroup of C_2 is 2.
order_sylow_C2 = 2

# The order of the Sylow 2-subgroup of G = A_5 x C_2 is the product
# of the orders of the Sylow 2-subgroups of A_5 and C_2.
sum_of_entries = order_sylow_A5 * order_sylow_C2

# The final equation is the product of the orders of these Sylow subgroups.
# We print each number involved in this final calculation.
print("The final sum is calculated from the equation:")
print(f"{order_sylow_A5} * {order_sylow_C2} = {sum_of_entries}")
print("\nThe sum of all the entries is:")
print(sum_of_entries)