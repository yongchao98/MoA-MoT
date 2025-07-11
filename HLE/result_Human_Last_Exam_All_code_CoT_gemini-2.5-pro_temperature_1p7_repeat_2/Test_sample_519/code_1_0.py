# Calculations for X1
# Dimension of the Hilbert scheme of points in affine space
degree_1 = 11
ambient_dim_1 = 3
dim1 = degree_1 * ambient_dim_1
profile1 = f"[S,s,{dim1}]"

# Calculations for X2
# Dimension of the quotient stack
ambient_dim_2 = 4
group_dim_2 = 1
dim2 = ambient_dim_2 - group_dim_2
profile2 = f"[DM,s,irr,{dim2}]"

# Calculations for X3
# Dimension of the Picard stack
genus_3 = 7
stabilizer_dim_3 = 1
dim3 = genus_3 - stabilizer_dim_3
profile3 = f"[A,s,{dim3}]"

# Combine the profiles into the final answer format
final_answer = f"{profile1} {profile2} {profile3}"
print(final_answer)