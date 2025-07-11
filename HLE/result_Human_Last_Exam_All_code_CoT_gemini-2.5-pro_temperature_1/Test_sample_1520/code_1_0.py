def get_generators(group_str):
    """
    Calculates the number of generators for a given group or a direct product of groups.
    Examples: "SU(3)", "SU(2) x U(1)"
    """
    # Handle direct product groups by summing the generators of their parts
    if "x" in group_str:
        parts = group_str.split("x")
        return sum(get_generators(part.strip()) for part in parts)
    
    group_str = group_str.strip().upper()
    
    # Calculate generators for SU(n)
    if group_str.startswith("SU("):
        try:
            n = int(group_str[3:-1])
            # The number of generators for SU(n) is n^2 - 1
            return n**2 - 1
        except (ValueError, IndexError):
            raise ValueError(f"Invalid SU(n) format: {group_str}")
    # Calculate generators for U(1)
    elif group_str == "U(1)":
        # The number of generators for U(1) is 1
        return 1
    else:
        raise ValueError(f"Unknown group type: {group_str}")

# Define the symmetry groups from the problem statement
initial_group = "SU(3)"
residual_group = "SU(2) x U(1)"

# --- Step-by-step Calculation ---

print("Analyzing the Spontaneous Symmetry Breaking: SU(3) -> SU(2) x U(1)")
print("=" * 70)

# Step 1: Calculate generators of the initial group G = SU(3)
dim_G = get_generators(initial_group)
n_G = 3
print(f"The initial group G is {initial_group}.")
print(f"The number of generators for SU({n_G}) is dim(G) = {n_G}^2 - 1 = {dim_G}.")
print("-" * 70)

# Step 2: Calculate generators of the residual group H = SU(2) x U(1)
dim_H = get_generators(residual_group)
dim_H_su2 = get_generators("SU(2)")
dim_H_u1 = get_generators("U(1)")
n_H_su2 = 2
print(f"The residual group H is {residual_group}.")
print(f"Its generators are the sum of dim(SU({n_H_su2})) and dim(U(1)).")
print(f"dim(H) = ({n_H_su2}^2 - 1) + 1 = {dim_H_su2} + {dim_H_u1} = {dim_H}.")
print("-" * 70)

# Step 3: Calculate the number of broken generators
num_broken_generators = dim_G - dim_H
print("The vacuum degeneracy is determined by the number of broken generators.")
print("Number of Broken Generators = dim(G) - dim(H)")
print(f"This gives us: {dim_G} - {dim_H} = {num_broken_generators} broken generators.")
print("-" * 70)

# Step 4: Conclusion
print("In a non-Abelian gauge theory, the Higgs mechanism causes one gauge boson to become massive for each broken generator.")
print(f"Therefore, a key consequence of this vacuum degeneracy is the emergence of {num_broken_generators} massive gauge bosons.")
print("\nBased on this calculation, the correct statement among the choices is 'Four massive gauge bosons'.")
