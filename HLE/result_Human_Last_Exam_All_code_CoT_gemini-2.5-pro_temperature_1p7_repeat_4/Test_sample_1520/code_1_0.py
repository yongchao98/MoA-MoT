import math

def get_su_generators(N):
    """Calculates the number of generators for SU(N)."""
    return N**2 - 1

def get_u1_generators():
    """Returns the number of generators for U(1)."""
    return 1

# 1. Define the initial and residual symmetry groups.
initial_group_name = "SU(3)"
residual_group_factors = ["SU(2)", "U(1)"]

# 2. Calculate generators for the initial group G = SU(3).
n_g = 3
generators_g = get_su_generators(n_g)

# 3. Calculate generators for the residual group H = SU(2) x U(1).
# These are the unbroken generators.
n_h1 = 2
generators_h1 = get_su_generators(n_h1)
generators_h2 = get_u1_generators()
generators_h_total = generators_h1 + generators_h2

# 4. Calculate the number of broken generators.
broken_generators = generators_g - generators_h_total

# 5. The number of massive gauge bosons is equal to the number of broken generators.
num_massive_bosons = broken_generators

# --- Output the results ---

print(f"Analysis of the symmetry breaking: {initial_group_name} -> {' x '.join(residual_group_factors)}")
print("-" * 60)
print(f"The number of generators for the initial group {initial_group_name} is: {n_g}^2 - 1 = {generators_g}")
print(f"The number of unbroken generators from the residual group SU(2) is: {n_h1}^2 - 1 = {generators_h1}")
print(f"The number of unbroken generators from the residual group U(1) is: {generators_h2}")
print(f"The total number of unbroken generators is the sum: {generators_h1} + {generators_h2} = {generators_h_total}")
print("-" * 60)
print("The number of broken generators is the difference between the initial and residual generators.")
print(f"Calculation: {generators_g} (total) - {generators_h_total} (unbroken) = {broken_generators} (broken)")
print()
print("According to the Higgs mechanism, each broken generator corresponds to a gauge boson that becomes massive.")
print(f"Therefore, the system has {num_massive_bosons} massive gauge bosons.")
print("-" * 60)
