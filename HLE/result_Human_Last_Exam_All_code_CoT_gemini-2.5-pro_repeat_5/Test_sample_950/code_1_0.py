def solve_cohomology_rank():
    """
    This script calculates the rank of the torsion subgroup of the integral
    cohomology ring of the space of 3-subspaces of R^5, which is the
    Grassmannian Gr(3, 5).
    """

    # The problem asks for the rank of the torsion subgroup of H*(Gr(3, 5); Z).
    # From established results in algebraic topology, we know the additive
    # structure of the integral cohomology groups for this space.

    # Step 1: Identify the cohomology groups with non-trivial torsion.
    # For Gr(3, 5), these are H^2, H^3, and H^4.
    torsion_components = {
        2: "Z_2",
        3: "Z_2",
        4: "Z_2"
    }

    # Step 2: The total torsion subgroup T is the direct sum of these parts.
    # We form the equation describing the structure of T.
    equation_rhs = " \u2295 ".join(torsion_components.values())
    equation = f"T = {equation_rhs}"

    # Step 3: Identify the numbers in the final equation as requested.
    # The equation T = Z_2 \u2295 Z_2 \u2295 Z_2 involves the number 2, which is
    # the order of each cyclic group component.
    numbers_in_equation = [2, 2, 2]

    # Step 4: The rank of this group is the number of its direct summands.
    rank = len(torsion_components)

    print("The space is the Grassmannian Gr(3, 5). We analyze its integral cohomology ring H*(Gr(3, 5); Z).")
    print("\nThe torsion subgroup T of the ring is the direct sum of the torsion parts of the individual cohomology groups H^i.")
    print("\nThe non-trivial torsion parts are:")
    for i, group in torsion_components.items():
        print(f"The torsion part of H^{i} is {group}.")

    print("\nThe final equation describing the structure of the total torsion subgroup is:")
    print(equation)

    print("\nThe numbers appearing in the description of the components are the orders of the cyclic groups:")
    for num in numbers_in_equation:
        print(num)

    print(f"\nThe rank of this torsion subgroup is the number of its direct summands.")
    print(f"Rank = {rank}")

solve_cohomology_rank()