def get_generators(group_name):
    """Calculates the number of generators for a given Lie group."""
    if group_name.startswith("SU"):
        N = int(group_name[2:])
        return N**2 - 1
    elif group_name == "U(1)":
        return 1
    else:
        raise ValueError("Unsupported group")

# Define the initial and residual groups
initial_group = "SU(3)"
residual_group_factors = ["SU(2)", "U(1)"]

# --- Step 1: Calculate generators of the initial group G = SU(3) ---
print("Spontaneous Symmetry Breaking: SU(3) -> SU(2) x U(1)\n")
print("Step 1: Calculate the number of generators for the initial group G = SU(3).")
dim_G = get_generators(initial_group)
print(f"The number of generators for SU(N) is N^2 - 1.")
print(f"For G = SU(3), the number of generators is 3^2 - 1 = {dim_G}.")
print("-" * 50)

# --- Step 2: Calculate generators of the residual group H = SU(2) x U(1) ---
print("Step 2: Calculate the number of generators for the residual group H = SU(2) x U(1).")
dim_H_su2 = get_generators(residual_group_factors[0])
dim_H_u1 = get_generators(residual_group_factors[1])
dim_H = dim_H_su2 + dim_H_u1
print(f"For the SU(2) factor, the number of generators is 2^2 - 1 = {dim_H_su2}.")
print(f"For the U(1) factor, the number of generators is {dim_H_u1}.")
print(f"The total number of unbroken generators for H is the sum: {dim_H_su2} + {dim_H_u1} = {dim_H}.")
print("-" * 50)

# --- Step 3: Calculate the number of broken generators ---
print("Step 3: Calculate the number of broken generators.")
broken_generators = dim_G - dim_H
print("The number of broken generators is dim(G) - dim(H).")
print(f"Number of broken generators = {dim_G} - {dim_H} = {broken_generators}.")
print("-" * 50)

# --- Step 4: Relate to physical consequences ---
print("Step 4: Determine the physical consequence of the broken generators.")
massive_gauge_bosons = broken_generators
print("In a non-Abelian gauge theory, the Higgs mechanism causes one gauge boson to become massive for each broken generator.")
print(f"Therefore, the number of massive gauge bosons is equal to the number of broken generators.")
print(f"Result: There are {massive_gauge_bosons} massive gauge bosons.")
print("\nThis corresponds to answer choice E.")
