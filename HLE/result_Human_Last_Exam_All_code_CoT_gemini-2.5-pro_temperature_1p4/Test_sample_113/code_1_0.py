# The number of codimension 2 boundary strata of M_bar(3,1) can be found
# by classifying the dual graphs of stable curves with g=3, n=1, and 2 nodes.
# The calculation is broken down by the number of components in the stable curve.

# Case 1: 1 component
# This corresponds to an irreducible curve of genus 1 with 2 nodes.
# The single marked point is on this component.
count_1_comp = 1
print(f"Number of strata with 1 component: {count_1_comp}")

# Case 2: 2 components
# The sum of component genera must be 2. There are two graph structures.
# Subcase 2.1: Cycle graph (two components connected at two points)
#   - Genera (2,0): Marked point must be on the g=0 component for stability.
count_2_comp_cycle_g20 = 1
#   - Genera (1,1): Components are symmetric.
count_2_comp_cycle_g11 = 1
# Subcase 2.2: Loop-and-bridge graph (one component has a node and is attached to the other)
#   - Genera (2,0): Loop must be on g=0 component. Marked point on either component gives a distinct stratum.
count_2_comp_loop_g20 = 2
#   - Genera (1,1): Marked point on either component (which are topologically distinct) gives a distinct stratum.
count_2_comp_loop_g11 = 2
count_2_comp = count_2_comp_cycle_g20 + count_2_comp_cycle_g11 + count_2_comp_loop_g20 + count_2_comp_loop_g11
print(f"Number of strata with 2 components: {count_2_comp_cycle_g20} + {count_2_comp_cycle_g11} + {count_2_comp_loop_g20} + {count_2_comp_loop_g11} = {count_2_comp}")


# Case 3: 3 components
# This is a chain of 3 components. The sum of their genera must be 3.
#   - Genera (2,1,0): The g=0 component must be in the middle with the marked point.
count_3_comp_g210 = 1
#   - Genera (1,1,1): Marked point can be on an end or the middle component.
count_3_comp_g111 = 2
count_3_comp = count_3_comp_g210 + count_3_comp_g111
print(f"Number of strata with 3 components: {count_3_comp_g210} + {count_3_comp_g111} = {count_3_comp}")

# Total number of strata
total_strata = count_1_comp + count_2_comp + count_3_comp
print(f"\nTotal number of codimension 2 boundary strata is the sum of these counts.")
print(f"Total = {count_1_comp} + {count_2_comp} + {count_3_comp} = {total_strata}")