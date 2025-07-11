import math

def solve_total_rank():
    """
    Calculates the total rank of the SO(4)-equivariant cohomology ring A
    in degrees up to 100.

    The Poincaré series of the ring A is P(t) = t^2 / (1 - t^4)^2.
    Expanding this gives P(t) = sum_{k=0 to inf} (k+1) * t^(4k+2).
    We need to sum the coefficients (k+1) for all non-negative integers k
    such that the degree 4k+2 is less than or equal to 100.
    """
    
    max_degree = 100
    total_rank = 0
    
    # The list `ranks_to_sum` will store the rank contribution from each allowed degree.
    ranks_to_sum = []

    k = 0
    while True:
        degree = 4 * k + 2
        if degree > max_degree:
            break
        
        # The rank in this degree is the coefficient (k+1)
        rank_in_degree = k + 1
        ranks_to_sum.append(rank_in_degree)
        total_rank += rank_in_degree
        k += 1

    # The final equation is the sum of all the individual ranks.
    # We construct a string to display this sum explicitly.
    equation_str = " + ".join(map(str, ranks_to_sum))
    
    print("The total rank is the sum of ranks in degrees of the form 4k+2 up to 100.")
    print("This corresponds to k from 0 to 24.")
    print("The final calculation is the sum of (k+1) for k in this range:")
    print(f"\n{equation_str} = {total_rank}")

solve_total_rank()