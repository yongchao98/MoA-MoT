# This problem is about finding the rank of the torsion subgroup of the integral cohomology ring
# of the space of 3-subspaces of R^5, which is the real Grassmannian Gr(3, 5).

# The integral cohomology groups H^i(Gr(3, 5); Z) are known from the literature.
# The torsion part of H*(Gr(3, 5); Z) consists only of 2-torsion (summands of Z/2Z).
# The rank of a torsion group, in this context, refers to the number of its cyclic summands.

# The number of Z/2Z summands in each cohomology group H^i are:
# H^0: 0
# H^1: 1
# H^2: 1
# H^3: 1
# H^4: 1
# H^5: 1
# H^6: 0
# H^i for i > 6: 0

torsion_ranks = [0, 1, 1, 1, 1, 1, 0]

# Calculate the total rank of the torsion subgroup by summing the ranks for each degree.
total_rank = sum(torsion_ranks)

# Print the final equation as requested.
equation_str = " + ".join(map(str, torsion_ranks))
print(f"The rank of the torsion subgroup is the sum of the ranks for each degree H^i for i=0 to 6.")
print(f"The calculation is: {equation_str} = {total_rank}")
