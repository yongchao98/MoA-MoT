import sys

def solve_task():
    """
    Calculates the maximal possible number of complex zeros for the given matrix system.
    The number depends on the size of the matrix, N.
    """
    
    # Prompt the user for the value of N
    try:
        n_str = input("Please enter the size of the matrix N (an integer >= 1): ")
        N = int(n_str)
        if N < 1:
            raise ValueError()
    except (ValueError, NameError):
        print("Invalid input. Please enter a positive integer for N.")
        return

    # For N=1, the only solution is real, so there are no complex zeros.
    if N == 1:
        max_zeros = 0
        explanation = f"For N = 1, the equation is A_11 + k_1 = 0, which has one real solution k_1 = -A_11. \nThus, the number of complex zeros (Re(k_j) != 0, Im(k_j) != 0) is 0."
    # For N>=2, the maximal number of complex zeros is N * 2^(N-1).
    else:
        # Using the formula N * 2^(N-1)
        max_zeros = N * (2**(N - 1))
        explanation = f"For N = {N}, the maximal number of complex zeros is given by the formula N * 2^(N-1).\nCalculation: {N} * 2^({N}-1) = {N} * {2**(N-1)} = {max_zeros}"

    print(explanation)
    print("\nFinal Answer:")
    print(f"The maximal possible number of complex zeros is {max_zeros}.")

solve_task()
