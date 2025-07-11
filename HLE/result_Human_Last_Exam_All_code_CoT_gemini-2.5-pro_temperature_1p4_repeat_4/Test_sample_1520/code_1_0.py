import math

def get_su_generators(N):
    """Calculates the number of generators for an SU(N) group."""
    return N**2 - 1

def get_u1_generators():
    """Returns the number of generators for a U(1) group."""
    return 1

# 1. Define initial and final groups
initial_group_n = 3
final_group_su_n = 2
final_group_u_n = 1

# 2. Calculate the number of generators for the initial group G = SU(3)
dim_G = get_su_generators(initial_group_n)

# 3. Calculate the number of generators for the residual group H = SU(2) x U(1)
dim_H_su2 = get_su_generators(final_group_su_n)
dim_H_u1 = get_u1_generators()
dim_H = dim_H_su2 + dim_H_u1

# 4. Calculate the number of broken generators
num_broken_generators = dim_G - dim_H

# 5. Print the step-by-step explanation and results
print("Analysis of the spontaneous symmetry breaking SU(3) -> SU(2) x U(1):")
print("-" * 70)
print(f"The number of generators for the initial group SU({initial_group_n}) is {initial_group_n}^2 - 1 = {dim_G}.")
print(f"The number of generators for the residual subgroup SU({final_group_su_n}) is {final_group_su_n}^2 - 1 = {dim_H_su2}.")
print(f"The number of generators for the residual subgroup U({final_group_u_n}) is {dim_H_u1}.")
print(f"The total number of generators for the residual symmetry H = SU(2) x U(1) is {dim_H_su2} + {dim_H_u1} = {dim_H}.")
print("\nIn a gauge theory, the number of massive gauge bosons equals the number of broken generators.")
print("The number of broken generators is the dimension of the initial group minus the dimension of the residual group.")
print("\nFinal Calculation:")
print(f"Number of Massive Gauge Bosons = (Generators of SU({initial_group_n})) - (Generators of SU({final_group_su_n}) x U({final_group_u_n}))")
# The final part of the code outputs each number in the final equation as requested
print(f"Result: {dim_G} - ({dim_H_su2} + {dim_H_u1}) = {num_broken_generators}")
print("\nThis result indicates there are 4 massive gauge bosons, which corresponds to option E.")
