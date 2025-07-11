def solve():
    """
    Calculates the largest value k such that for any valid arrangement of k diamonds
    on a 2024x2024 grid, a move is always possible.

    The logic is based on finding the size of the smallest "stuck" arrangement, m(N),
    and then calculating k = m(N) - 1.

    For an N x N grid where N is a multiple of 4, the minimum size of a stuck
    arrangement is conjectured to be m(N) = (N/2)^2 + N - 1.
    Therefore, the desired value is k = (N/2)^2 + N - 2.
    """
    N = 2024
    
    # The size of the smallest stuck arrangement m(N) is (N/2)^2 + N - 1
    # The required value k is m(N) - 1
    k = (N // 2)**2 + N - 2
    
    # Output the final equation and the result
    L = N // 2
    print(f"The formula for the largest value k is (N/2)^2 + N - 2.")
    print(f"For N = {N}:")
    print(f"k = ({N}//2)^2 + {N} - 2")
    print(f"k = {L}^2 + {N - 2}")
    print(f"k = {L**2} + {N - 2}")
    print(f"k = {k}")

solve()