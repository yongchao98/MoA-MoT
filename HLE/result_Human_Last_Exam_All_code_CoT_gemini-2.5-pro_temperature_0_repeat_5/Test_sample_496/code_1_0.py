def solve():
    """
    Calculates the total rank of the equivariant cohomology ring A
    for degrees up to 100.
    """
    limit = 100
    total_rank = 0

    # Part 1: Ranks for degrees of the form 4k
    k_max_1 = limit // 4
    ranks_1 = list(range(1, k_max_1 + 1))
    sum_1 = sum(ranks_1)
    
    # Part 2: Ranks for degrees of the form 4k+2
    k_max_2 = (limit - 2) // 4
    ranks_2 = list(range(1, k_max_2 + 1))
    sum_2 = sum(ranks_2)

    total_rank = sum_1 + sum_2

    # Output the explanation and the final equation
    print("The total rank is the sum of ranks for degrees of the form 4k and 4k+2 up to 100.")
    
    # Print the first part of the equation
    print("Sum of ranks for degrees of the form 4k (k=1 to 25):")
    print(" + ".join(map(str, ranks_1)) + f" = {sum_1}")

    # Print the second part of the equation
    print("\nSum of ranks for degrees of the form 4k+2 (k=1 to 24):")
    print(" + ".join(map(str, ranks_2)) + f" = {sum_2}")

    # Print the final calculation
    print(f"\nTotal Rank = {sum_1} + {sum_2} = {total_rank}")

solve()