def solve_total_rank():
    """
    Calculates the total rank of the SO(4)-equivariant cohomology ring of SO(4)\X
    in degrees up to 100.
    """
    
    # The degrees for which the rank is non-zero are of the form 4k.
    # For a degree to be <= 100, we have 4k <= 100, so k <= 25.
    # We sum k from 0 to 25.
    limit_k = 25
    
    # The rank in degree 4k is 2*(k+1).
    ranks = [2 * (k + 1) for k in range(limit_k + 1)]
    
    total_rank = sum(ranks)
    
    # Building the equation string showing the sum.
    equation_str = " + ".join(map(str, ranks))
    
    print("The total rank is the sum of ranks in non-zero cohomology groups up to degree 100.")
    print("Non-zero ranks are in degrees 4k, where k=0, 1, ..., 25.")
    print(f"The final sum is: {equation_str}")
    print(f"Total Rank = {total_rank}")

solve_total_rank()