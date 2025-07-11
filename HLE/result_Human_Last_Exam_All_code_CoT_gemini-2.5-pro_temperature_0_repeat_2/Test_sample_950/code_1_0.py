# Step 1: Define the ranks of the torsion subgroups for each cohomology group H^i.
# Based on the known results for the integral cohomology of Gr_3(R^5).
# H^0: Z (rank 0)
# H^1: 0 (rank 0)
# H^2: Z/2Z (rank 1)
# H^3: Z/2Z (rank 1)
# H^4: Z + Z/2Z (rank 1)
# H^5: Z/2Z + Z/2Z (rank 2)
# H^6: Z/2Z (rank 1)
# H^i for i > 6 are 0.

torsion_ranks = {
    0: 0,
    1: 0,
    2: 1,
    3: 1,
    4: 1,
    5: 2,
    6: 1
}

# Step 2: Calculate the total rank by summing the ranks for each degree.
total_rank = sum(torsion_ranks.values())

# Step 3: Print the calculation step-by-step and the final result.
ranks_list = [str(r) for r in torsion_ranks.values() if r > 0]
equation = " + ".join(ranks_list)
print(f"The rank of the torsion subgroup is the sum of the ranks of the torsion part of each cohomology group H^i.")
print(f"The ranks for i=2, 3, 4, 5, 6 are 1, 1, 1, 2, 1 respectively.")
print(f"The total rank is {equation} = {total_rank}")
