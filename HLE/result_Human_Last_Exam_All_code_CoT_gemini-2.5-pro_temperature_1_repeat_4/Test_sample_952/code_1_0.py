import math

def solve():
    """
    Solves the diamond arrangement puzzle.

    The problem asks for the largest value K such that for every valid arrangement
    of K diamonds, at least one diamond can be moved to an adjacent cell.
    This is equivalent to finding M - 1, where M is the size of the smallest
    "stuck" arrangement (where no diamond can be moved).

    A "stuck" arrangement is a specific type of "maximal independent set". The size
    of the smallest maximal independent set on an N x N grid is given by the
    formula ceil(N*N / 5) for large N. We assume that this is the size M we
    are looking for.

    The grid size is N = 2024.
    """
    N = 2024
    
    # Calculate the total number of cells
    total_cells = N * N
    
    # Calculate the size of the smallest maximal independent set, M
    # M = ceil(N^2 / 5)
    M = math.ceil(total_cells / 5)
    
    # The largest value K for which any arrangement is not "stuck" is M - 1.
    K = M - 1
    
    print(f"The side length of the square table is N = {N}.")
    print(f"The total number of cells is N * N = {N} * {N} = {total_cells}.")
    print(f"The size of the smallest 'stuck' (minimal maximal) arrangement is M = ceil({total_cells} / 5) = {M}.")
    print(f"The largest value K such that any arrangement of K diamonds has a possible move is K = M - 1.")
    print(f"K = {M} - 1 = {K}.")
    print("\nFinal Answer:")
    print(K)

solve()