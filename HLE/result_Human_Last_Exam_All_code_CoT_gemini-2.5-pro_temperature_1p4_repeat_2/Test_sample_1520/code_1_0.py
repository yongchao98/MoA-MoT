import math

def calculate_generators(group_name):
    """Calculates the number of generators for SU(N) or U(1) groups."""
    if group_name.startswith("SU"):
        n = int(group_name[2])
        return n**2 - 1
    elif group_name == "U(1)":
        return 1
    else:
        return 0

# Define the groups
initial_group = "SU(3)"
residual_group_su2 = "SU(2)"
residual_group_u1 = "U(1)"

# Calculate generators for the initial group G = SU(3)
g_initial = calculate_generators(initial_group)
n_initial = int(initial_group[2])

# Calculate generators for the residual group H = SU(2) x U(1)
g_residual_su2 = calculate_generators(residual_group_su2)
n_residual_su2 = int(residual_group_su2[2])
g_residual_u1 = calculate_generators(residual_group_u1)
g_residual_total = g_residual_su2 + g_residual_u1

# The number of broken generators equals the number of massive gauge bosons
num_massive_bosons = g_initial - g_residual_total

# Print the detailed calculation
print("In a non-Abelian gauge theory, the number of massive gauge bosons is determined by the number of broken generators.")
print("This is calculated as the difference between the generators of the initial group G and the residual group H.")
print("\nCalculation:")
print(f"Number of massive gauge bosons = (Generators of G) - (Generators of H)")
print(f"                               = (Generators of {initial_group}) - (Generators of {residual_group_su2} + Generators of {residual_group_u1})")
print(f"                               = ({n_initial}^2 - 1) - (({n_residual_su2}^2 - 1) + 1)")
print(f"                               = ({int(math.pow(n_initial, 2))} - 1) - (({int(math.pow(n_residual_su2, 2))} - 1) + {g_residual_u1})")
print(f"                               = {g_initial} - ({g_residual_su2} + {g_residual_u1})")
print(f"                               = {g_initial} - {g_residual_total}")
print(f"                               = {num_massive_bosons}")
print(f"\nTherefore, the symmetry breaking results in {num_massive_bosons} massive gauge bosons.")
