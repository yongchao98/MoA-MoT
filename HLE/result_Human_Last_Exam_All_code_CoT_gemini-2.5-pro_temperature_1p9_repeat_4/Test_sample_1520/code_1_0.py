def dim_su_n(n):
    """Calculates the number of generators for the SU(N) group."""
    return n**2 - 1

def dim_u_1():
    """Returns the number of generators for the U(1) group."""
    return 1

# 1. Calculate the number of generators for the initial group G = SU(3)
initial_generators = dim_su_n(3)
print(f"The initial symmetry group SU(3) has {initial_generators} generators.")

# 2. Calculate the number of generators for the residual group H = SU(2) x U(1)
su2_generators = dim_su_n(2)
u1_generators = dim_u_1()
residual_generators = su2_generators + u1_generators
print(f"The residual symmetry group SU(2) x U(1) has {su2_generators} + {u1_generators} = {residual_generators} generators.")

# 3. Calculate the number of broken generators
broken_generators = initial_generators - residual_generators
print("\nThe number of broken generators is the difference between the initial and residual generators.")
print(f"Final Equation: {initial_generators} - {residual_generators} = {broken_generators}")

print("\nThis result corresponds to the number of broken generators, which in a gauge theory, also equals the number of massive gauge bosons that acquire mass via the Higgs mechanism.")
print("Each broken generator corresponds to a direction in the space of vacuum states, defining the vacuum degeneracy.")
