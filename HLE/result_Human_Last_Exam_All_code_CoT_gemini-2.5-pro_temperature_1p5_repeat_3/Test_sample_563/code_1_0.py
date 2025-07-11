# The number of isomorphism classes of automorphism groups for Riemann surfaces
# of a given genus is a known result from algebraic geometry.
# These values are based on catalogs published in mathematical research papers.

# For genus g=2, the number of groups is 12.
# All genus 2 curves are hyperelliptic, and their automorphism groups are well-classified.
num_groups_g2 = 12

# For genus g=3, the number of groups is 36.
# This involves classifying groups for both hyperelliptic and non-hyperelliptic curves.
num_groups_g3 = 36

# For genus g=4, the number of groups is 23.
# This classification is the most complex of the three and is detailed in
# research by M. Conder and others.
num_groups_g4 = 23

# The final result is a list containing these three numbers.
# We print each number within the final list structure as requested.
print(f"[{num_groups_g2},{num_groups_g3},{num_groups_g4}]")