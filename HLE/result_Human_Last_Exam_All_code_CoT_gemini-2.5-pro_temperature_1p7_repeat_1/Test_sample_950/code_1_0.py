# Step 1: Define the ranks of the torsion subgroups in each degree.
# These ranks are calculated by finding the rank of the linear map
# Sq^1: H^i(G_2(R^5); Z_2) -> H^{i+1}(G_2(R^5); Z_2) for i = 0 to 6.
# Torsion rank in H^{i+1} is the rank of Sq^1 on H^i.

# Torsion rank in H^0 is zero by definition.
rank_T0 = 0
# Torsion rank in H^1 is zero because H^0 has no torsion. This corresponds to the rank of Sq^1 on H^{-1}, which is 0.
rank_T1 = 0
# Torsion rank in H^2 is rank(Sq^1: H^1 -> H^2)
rank_T2 = 1
# Torsion rank in H^3 is rank(Sq^1: H^2 -> H^3)
rank_T3 = 1
# Torsion rank in H^4 is rank(Sq^1: H^3 -> H^4)
rank_T4 = 2
# Torsion rank in H^5 is rank(Sq^1: H^4 -> H^5)
rank_T5 = 0
# Torsion rank in H^6 is rank(Sq^1: H^5 -> H^6)
rank_T6 = 1
# Torsion rank for degrees > 6 is 0 as the cohomology groups are zero.
ranks = [rank_T0, rank_T1, rank_T2, rank_T3, rank_T4, rank_T5, rank_T6]

# Step 2: Sum the ranks to find the total rank of the torsion subgroup.
total_rank = sum(ranks)

# Step 3: Print the calculation step by step.
print("The total rank of the torsion subgroup is the sum of the ranks of the torsion subgroups in each cohomology degree.")
print("The rank of the torsion in degree i, rank(T^i), is calculated.")
print(f"rank(T^0) = {ranks[0]}")
print(f"rank(T^1) = {ranks[1]}")
print(f"rank(T^2) = {ranks[2]}")
print(f"rank(T^3) = {ranks[3]}")
print(f"rank(T^4) = {ranks[4]}")
print(f"rank(T^5) = {ranks[5]}")
print(f"rank(T^6) = {ranks[6]}")

# Create the equation string
equation_parts = [str(r) for r in ranks if r > 0]
equation_str = " + ".join(equation_parts)

print(f"\nTotal rank = {equation_str} = {total_rank}")
