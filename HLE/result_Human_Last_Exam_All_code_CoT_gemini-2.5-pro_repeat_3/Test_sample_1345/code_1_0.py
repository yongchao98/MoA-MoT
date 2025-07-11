def solve():
    """
    Calculates the maximal possible number of complex zeros based on the analysis above.
    The formula derived is M(N) = 0 for N=1, and M(N) = N * 2^(N-1) for N > 1.
    This script prints the result of this formula for a few example values of N.
    """
    print("The maximal possible number of complex zeros is given by the formula:")
    print("M(N) = 0, if N = 1")
    print("M(N) = N * 2^(N-1), if N > 1")
    print("\nThis can be demonstrated for a few values of N:")
    
    for n_val in range(1, 6):
        if n_val == 1:
            result = 0
            print(f"For N = {n_val}, the maximal number is {result}")
        else:
            # The final equation is N * 2^(N-1)
            power = n_val - 1
            term2 = 2**power
            result = n_val * term2
            
            # Outputting each number in the final equation as requested
            print(f"For N = {n_val}, the maximal number is {n_val} * 2^({power}) = {result}")

solve()