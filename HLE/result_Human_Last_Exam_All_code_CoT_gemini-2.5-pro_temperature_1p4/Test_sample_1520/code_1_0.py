import math

def get_su_n_generators(n):
  """Calculates the number of generators for an SU(N) group."""
  if n < 2:
    return 0
  return n**2 - 1

def get_u_1_generators():
  """Returns the number of generators for a U(1) group."""
  return 1

# 1. Define the initial and residual groups from the problem statement.
initial_group_n = 3
residual_group_su_n = 2

# 2. Calculate the number of generators for the initial group G = SU(3).
generators_g = get_su_n_generators(initial_group_n)

# 3. Calculate the number of generators for the residual group H = SU(2) x U(1).
generators_h = get_su_n_generators(residual_group_su_n) + get_u_1_generators()

# 4. Calculate the number of broken generators, which defines the vacuum degeneracy.
broken_generators = generators_g - generators_h

# In a gauge theory, this is equal to the number of massive gauge bosons.
massive_gauge_bosons = broken_generators

print(f"The spontaneous symmetry breaking is given by SU({initial_group_n}) -> SU({residual_group_su_n}) x U(1).\n")
print(f"Number of generators in the initial group SU({initial_group_n}): {generators_g}")
print(f"Number of generators in the residual group SU({residual_group_su_n}) x U(1): {generators_h}")
print("-" * 40)
print("The number of massive gauge bosons is equal to the number of broken generators.")
print("This is calculated as the difference between the generators of the initial and residual groups.")
print("\nFinal Equation:")
print(f"Number of massive gauge bosons = {generators_g} (from SU({initial_group_n})) - {generators_h} (from SU({residual_group_su_n})xU(1)) = {massive_gauge_bosons}")
print("\nThis result corresponds to the statement: Four massive gauge bosons.")