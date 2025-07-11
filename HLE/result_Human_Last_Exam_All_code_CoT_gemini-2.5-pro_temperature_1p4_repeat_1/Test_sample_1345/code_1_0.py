import sys

def solve():
    """
    Calculates the maximal possible number of complex zeros for the given matrix B(k).
    The user is prompted to enter the size of the matrix, N.
    """
    try:
        # Prompt the user to enter the value of N
        N_str = input("Enter the dimension of the matrix N: ")
        N = int(N_str)
        if N < 1:
            print("Error: N must be a positive integer.", file=sys.stderr)
            return

    except ValueError:
        print(f"Error: Invalid input '{N_str}'. Please enter an integer.", file=sys.stderr)
        return

    # Case N=1
    # The equation is A11 + k1 = 0. Since A11 is real, k1 must be real.
    # So, there are no complex zeros.
    if N == 1:
        result = 0
        print("\nFor N=1, the equation is A_11 + k_1 = 0.")
        print("Since A_11 is real, k_1 is always real.")
        print("Maximum number of complex zeros = 0")

    # Case N >= 2
    # The number of solutions is given by the degree of the final polynomial in k1,
    # which is N * 2^(N-1).
    else:
        # Calculate the result using the formula
        power = N - 1
        term = 2**power
        result = N * term
        
        # Print the final equation with all the numbers
        print(f"\nFor N >= 2, the maximal number of complex zeros is given by the formula N * 2^(N-1).")
        print(f"Calculation for N = {N}:")
        print(f"{N} * 2^({N} - 1) = {N} * 2^{power} = {N} * {term} = {result}")
        print(f"Maximum number of complex zeros = {result}")

if __name__ == "__main__":
    solve()