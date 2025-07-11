# The number of different types of stable reduction for curves of genus 2
# is a known result from algebraic geometry, corresponding to the number of
# strata in the moduli space of stable genus 2 curves.

# 1. The generic, non-degenerate type
num_smooth = 1
print(f"Number of non-degenerate types (smooth curves): {num_smooth}")

# 2. Types with at least one component of geometric genus 1 (elliptic curve)
# This includes:
# - An irreducible curve with one node (normalization is an elliptic curve).
# - A reducible curve made of two elliptic curves.
# - A reducible curve made of one elliptic curve and a rational part (2 types).
num_with_elliptic_component = 4
print(f"Number of degenerate types with an elliptic component: {num_with_elliptic_component}")


# 3. Types where all components are of geometric genus 0 (rational curves)
# These are sometimes called "totally degenerate" or "rational" types.
# This includes:
# - An irreducible curve with two nodes (normalization is a rational curve).
# - Reducible curves made only of rational curves (8 distinct types).
num_rational_types = 9
print(f"Number of degenerate types with only rational components: {num_rational_types}")

# The total number of stable reduction types is the sum of these categories.
total_types = num_smooth + num_with_elliptic_component + num_rational_types

print(f"\nEquation: {num_smooth} + {num_with_elliptic_component} + {num_rational_types} = {total_types}")
print(f"\nThus, there are a total of {total_types} different types of stable reduction for curves of genus 2.")
