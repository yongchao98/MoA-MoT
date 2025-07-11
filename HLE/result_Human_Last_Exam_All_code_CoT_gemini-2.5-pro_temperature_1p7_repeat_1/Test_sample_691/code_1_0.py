# In group theory, Z represents the group of integers.
# The fundamental group of a torus (T^2) is the direct product Z x Z.
torus_group = "(Z x Z)"

# The final space is topologically equivalent to the wedge sum of two tori.
# The fundamental group of a wedge sum is the free product of the individual groups.
# The free product operator is denoted by '*'.
final_group_expression = f"{torus_group} * {torus_group}"

# Print the final symbolic representation of the fundamental group.
print(final_group_expression)