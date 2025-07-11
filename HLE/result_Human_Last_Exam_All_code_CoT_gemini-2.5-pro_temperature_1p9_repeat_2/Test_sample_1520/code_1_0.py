import re

def get_generators(group_str: str) -> int:
    """
    Calculates the number of generators for a Lie group specified by a string.
    Supports SU(N) and U(1) groups.
    """
    # Check for SU(N) pattern
    if match := re.match(r"SU\((\d+)\)", group_str, re.IGNORECASE):
        n = int(match.group(1))
        return n**2 - 1
    # Check for U(1) pattern
    elif re.match(r"U\(1\)", group_str, re.IGNORECASE):
        return 1
    else:
        # Raise an error for unsupported group types
        raise ValueError(f"Unsupported group format: {group_str}")

def analyze_symmetry_breaking():
    """
    Analyzes the spontaneous symmetry breaking SU(3) -> SU(2) x U(1)
    and prints the step-by-step reasoning.
    """
    # Define the symmetry breaking scenario
    initial_group = "SU(3)"
    residual_group_parts = ["SU(2)", "U(1)"]

    print("Analyzing the spontaneous symmetry breaking SU(3) -> SU(2) x U(1):")
    print("-" * 60)

    # Step 1: Calculate generators for the initial group G = SU(3)
    dim_G = get_generators(initial_group)
    print(f"1. The initial symmetry group is G = {initial_group}.")
    print(f"   The number of total generators for {initial_group} (N=3) is 3^2 - 1 = {dim_G}.")
    print("   This corresponds to the total number of gauge bosons before breaking.")
    print("\n")

    # Step 2: Calculate generators for the residual group H = SU(2) x U(1)
    dim_H_su2 = get_generators(residual_group_parts[0])
    dim_H_u1 = get_generators(residual_group_parts[1])
    dim_H_total = dim_H_su2 + dim_H_u1
    print(f"2. The residual (unbroken) symmetry group is H = {residual_group_parts[0]} x {residual_group_parts[1]}.")
    print(f"   The {residual_group_parts[0]} part has 2^2 - 1 = {dim_H_su2} unbroken generators.")
    print(f"   The {residual_group_parts[1]} part has {dim_H_u1} unbroken generator.")
    print(f"   The total number of unbroken generators is {dim_H_su2} + {dim_H_u1} = {dim_H_total}.")
    print("   These correspond to the gauge bosons that remain massless.")
    print("\n")

    # Step 3: Calculate the number of broken generators
    num_broken = dim_G - dim_H_total
    print("3. The 'vacuum degeneracy' arises from the symmetries that are broken.")
    print(f"   The number of broken generators = (Total generators of G) - (Total generators of H).")
    print(f"   This gives us {dim_G} - {dim_H_total} = {num_broken} broken generators.")
    print("\n")

    # Step 4: Relate broken generators to massive gauge bosons via the Higgs mechanism
    print("4. In a non-Abelian gauge theory, the Higgs mechanism takes effect:")
    print("   Each broken generator gives rise to a massive gauge boson.")
    print(f"   Therefore, {num_broken} broken generators imply the existence of {num_broken} massive gauge bosons.")
    print("-" * 60)
    print("This result, 'Four massive gauge bosons', uniquely characterizes a key physical")
    print("consequence of the vacuum degeneracy in this specific gauge theory.")
    print("\nFinal equation showing the calculation:")
    print(f"{dim_G} - ({dim_H_su2} + {dim_H_u1}) = {num_broken}")

if __name__ == '__main__':
    analyze_symmetry_breaking()