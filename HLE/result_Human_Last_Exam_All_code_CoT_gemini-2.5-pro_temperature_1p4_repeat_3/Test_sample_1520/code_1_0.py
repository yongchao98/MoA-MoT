import math

def get_su_generators(N):
  """Calculates the number of generators for the SU(N) group."""
  return N**2 - 1

def get_u1_generators():
  """Returns the number of generators for the U(1) group."""
  return 1

# 1. Define the initial and residual groups
initial_group_n = 3
residual_group_su_n = 2

# 2. Calculate generators for the initial group G = SU(3)
generators_G = get_su_generators(initial_group_n)
print(f"The initial symmetry group is SU({initial_group_n}).")
print(f"Number of generators for SU({initial_group_n}) = {initial_group_n}**2 - 1 = {generators_G}")
print("-" * 30)

# 3. Calculate generators for the residual group H = SU(2) x U(1)
generators_H_su2 = get_su_generators(residual_group_su_n)
generators_H_u1 = get_u1_generators()
generators_H_total = generators_H_su2 + generators_H_u1
print(f"The residual symmetry group is SU({residual_group_su_n}) x U(1).")
print(f"Number of generators for SU({residual_group_su_n}) = {residual_group_su_n}**2 - 1 = {generators_H_su2}")
print(f"Number of generators for U(1) = {generators_H_u1}")
print(f"Total number of unbroken generators = {generators_H_su2} + {generators_H_u1} = {generators_H_total}")
print("-" * 30)

# 4. Calculate the number of broken generators
broken_generators = generators_G - generators_H_total
print("The number of broken generators is the difference between the initial and residual generators.")
print(f"Number of broken generators = {generators_G} - {generators_H_total} = {broken_generators}")
print("-" * 30)

# 5. Relate to massive gauge bosons
massive_gauge_bosons = broken_generators
print("In a non-Abelian gauge theory, the number of massive gauge bosons is equal to the number of broken generators.")
print(f"Therefore, the number of massive gauge bosons is {massive_gauge_bosons}.")
