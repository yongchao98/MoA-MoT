# This script outlines the properties of a 1D Fermi-Hubbard system with
# two-body losses in the infinite-time limit.

# Based on the analysis, the final state is the vacuum, which is characterized
# by zero particles, zero tunneling, and zero losses.
final_state_properties = [1, 2, 3]

# The dynamics leading to this state favor the formation of anti-ferromagnetic
# correlations to minimize the probability of creating doubly-occupied sites
# that lead to loss.
transient_dynamics_property = 5

# The full set of properties is the union of these two sets.
all_property_numbers = sorted(final_state_properties + [transient_dynamics_property])

print("The analysis shows that the system is defined by properties related to its final vacuum state and the transient dynamics.")
print("\nThe final set of identifying properties corresponds to the following numbers:")

# This loop prints each number of the final property set, as requested.
for number in all_property_numbers:
    print(number)
