# Dimension calculation and profile generation for X_1
dim_X1_n = 3
dim_X1_d = 11
dim_X1 = dim_X1_n * dim_X1_d
# Properties: [Scheme, separated, dimension]
# It is reducible, so 'irr' is not included. It is not universally closed.
prop_X1 = f"[S, s, {dim_X1}]"

# Dimension calculation and profile generation for X_2
dim_X2_base = 4
dim_X2_group = 1
dim_X2 = dim_X2_base - dim_X2_group
# Properties: [Algebraic Space, separated, irreducible, dimension]
# It is not universally closed.
prop_X2 = f"[S, s, irr, {dim_X2}]"

# Dimension calculation and profile generation for X_3
dim_X3_g = 7
dim_X3 = dim_X3_g
# Properties: [Algebraic stack, separated, dimension]
# It is reducible and not universally closed.
prop_X3 = f"[A, s, {dim_X3}]"

# The problem asks to output each number in the final equation.
# Here we print the equations for the dimensions.
print(f"Dimension of X_1: {dim_X1_n} * {dim_X1_d} = {dim_X1}")
print(f"Dimension of X_2: {dim_X2_base} - {dim_X2_group} = {dim_X2}")
print(f"Dimension of X_3: genus g = {dim_X3}")

# Combine the profiles into the final answer format
final_answer = f"{prop_X1} {prop_X2} {prop_X3}"
print("\nFinal Profiles:")
print(final_answer)