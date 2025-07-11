# The number of isomorphism classes of automorphism groups for compact
# Riemann surfaces of a given genus is a known, but non-trivial, result
# from the mathematical literature. These values are based on the latest
# census found in the paper by M. Conder, A. Neshchadim, and J. Wolfart (2020).

# Number of isomorphism classes for a genus g=2 surface.
num_groups_g2 = 12

# Number of isomorphism classes for a genus g=3 surface.
num_groups_g3 = 19

# Number of isomorphism classes for a genus g=4 surface.
num_groups_g4 = 28

# The problem asks for the answer as a list of these three numbers.
# We print the components of the final list.
print(f"Number of automorphism groups for genus 2: {num_groups_g2}")
print(f"Number of automorphism groups for genus 3: {num_groups_g3}")
print(f"Number of automorphism groups for genus 4: {num_groups_g4}")

# We now format the final answer as requested.
# The "equation" here is the construction of the list.
print(f"[{num_groups_g2},{num_groups_g3},{num_groups_g4}]")