def solve():
    """
    This function finds all integers n between 7 and 55 (inclusive) for which
    it is possible to leave a single gift on an n x n x n x n x n hypercube.

    The condition for this to be possible is that n must be congruent to
    1 or 13 modulo 14.
    """
    
    valid_n = []
    for n in range(7, 56):
        # Check if n = 14k + 1 or n = 14k + 13
        if n % 14 == 1 or n % 14 == 13:
            valid_n.append(n)
            
    # The problem asks to output each number. We print them separated by commas.
    # The list is already in increasing order.
    print(*valid_n, sep=', ')

solve()