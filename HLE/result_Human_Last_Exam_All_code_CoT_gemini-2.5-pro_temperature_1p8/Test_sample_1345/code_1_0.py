import sys

def solve():
    """
    Calculates the maximal possible number of complex zeros for the given matrix problem.
    The user is prompted to enter the size of the matrix, N.
    """
    try:
        # Take N as input from the user
        N_str = input("Enter the size of the matrix N (an integer >= 2): ")
        N = int(N_str)
        if N < 2:
            print("Error: N must be an integer greater than or equal to 2.", file=sys.stderr)
            return
    except ValueError:
        print("Error: Invalid input. Please enter an integer.", file=sys.stderr)
        return

    # Calculate the maximal number of complex zeros based on whether N is even or odd
    if N % 2 == 0:
        # If N is even
        max_zeros = N * (2**(N - 1))
        # The equation for N is N * 2^(N-1)
        equation_str = f"{N} * 2^({N}-1)"
    else:
        # If N is odd
        max_zeros = (N - 1) * (2**(N - 1))
        # The equation for N is (N-1) * 2^(N-1)
        equation_str = f"({N}-1) * 2^({N}-1)"

    # Print the explanation and the final result
    print("The maximal possible number of complex zeros is given by the formula:")
    print(f"({N} - ({N} % 2)) * 2^({N}-1)")
    print("\nFor N = {}, the calculation is:".format(N))
    print(f"{equation_str} = {max_zeros}")
    print("\nFinal Answer:")
    print(max_zeros)

if __name__ == "__main__":
    solve()
