import math

def calculate_symmetry_breaking(initial_group_n, residual_su_n, residual_u1):
    """
    Calculates the properties of a spontaneous symmetry breaking G -> H.

    Args:
        initial_group_n (int): The N for the initial SU(N) group.
        residual_su_n (int): The N for the residual SU(N) part of the group.
        residual_u1 (bool): True if there is a U(1) factor in the residual group.
    """

    # 1. Calculate generators for the initial group G = SU(initial_group_n)
    dim_g = initial_group_n**2 - 1
    print(f"The initial symmetry group is G = SU({initial_group_n}).")
    print(f"The number of generators for G is {initial_group_n}^2 - 1 = {dim_g}.")
    print("-" * 50)

    # 2. Calculate generators for the residual group H = SU(residual_su_n) x U(1)
    dim_h_su = residual_su_n**2 - 1
    dim_h_u1 = 1 if residual_u1 else 0
    dim_h = dim_h_su + dim_h_u1
    print(f"The residual symmetry group is H = SU({residual_su_n}) x U(1).")
    print(f"The number of unbroken generators is:")
    print(f"dim(SU({residual_su_n})) + dim(U(1)) = ({residual_su_n}^2 - 1) + {dim_h_u1} = {dim_h_su} + {dim_h_u1} = {dim_h}.")
    print("-" * 50)

    # 3. Calculate the number of broken generators
    num_broken = dim_g - dim_h
    print("The vacuum degeneracy is determined by the number of broken generators.")
    print(f"Number of broken generators = dim(G) - dim(H) = {dim_g} - {dim_h} = {num_broken}.")
    print("-" * 50)

    # 4. Relate broken generators to massive gauge bosons via Higgs mechanism
    num_massive_bosons = num_broken
    print("In a gauge theory, the Higgs mechanism dictates that the number of massive gauge")
    print("bosons is equal to the number of broken generators.")
    print(f"Therefore, the physical consequence of this symmetry breaking is the emergence of {num_massive_bosons} massive gauge bosons.")
    print(f"Final Equation: ({initial_group_n}^2 - 1) - (({residual_su_n}^2 - 1) + 1) = {dim_g} - ({dim_h_su} + 1) = {num_broken} massive gauge bosons.")

# For the specific case SU(3) -> SU(2) x U(1)
calculate_symmetry_breaking(initial_group_n=3, residual_su_n=2, residual_u1=True)

<<<E>>>