def solve_tori_count():
    """
    Calculates the number of Fq-rational maximal tori for a reductive group
    of type E8 over a finite field Fq with q elements.
    """

    # Group properties for type E8
    group_type = "E8"
    dimension = 248
    rank = 8

    print(f"Finding the number of rational maximal tori for a group of type {group_type} over the finite field Fq.")
    print("The number is given by the formula q^N, where N is the number of positive roots.")
    print("-" * 50)

    # Step 1: Calculate the total number of roots
    print("Step 1: Calculate the total number of roots for the E8 root system.")
    print(f"The dimension of the Lie algebra of type {group_type} is {dimension}.")
    print(f"The rank of the Lie algebra of type {group_type} is {rank}.")
    total_roots = dimension - rank
    print(f"The total number of roots = Dimension - Rank = {dimension} - {rank} = {total_roots}.")
    print("-" * 50)

    # Step 2: Calculate the number of positive roots (N)
    print("Step 2: Calculate the number of positive roots, N.")
    positive_roots = total_roots // 2
    print(f"The number of positive roots (N) is half the total number of roots.")
    print(f"N = Total Roots / 2 = {total_roots} / 2 = {positive_roots}.")
    print("-" * 50)

    # Step 3: State the final answer
    print("Step 3: State the final formula for the number of tori.")
    print(f"The number of Fq-rational maximal tori is q^N.")
    print(f"Substituting N = {positive_roots}, we get the final answer.")
    print(f"\nFinal Answer Equation: Number of Tori = q^{positive_roots}")

solve_tori_count()
