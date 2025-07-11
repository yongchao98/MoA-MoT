import sys

def solve():
    """
    Calculates the maximal possible number of complex zeros for the determinant
    of the given N x N matrix B(k).
    """
    try:
        # Prompting the user to enter the size of the matrix, N.
        n_str = input("Enter the size of the matrix N (an integer >= 1): ")
        n = int(n_str)
        if n < 1:
            print("Error: N must be a positive integer.", file=sys.stderr)
            return

        # For N=1, the only zero is real, so there are no complex zeros.
        # For N > 1, the maximal number of complex zeros is given by 2*N*(N-1).
        if n == 1:
            result = 0
            print(f"For N = 1, the maximal number of complex zeros is 0.")
            print(f"The equation is: 2 * 1 * (1 - 1) = {result}")

        else:
            result = 2 * n * (n - 1)
            print(f"For N = {n}, the maximal possible number of complex zeros is {result}.")
            print(f"The equation is: 2 * {n} * ({n} - 1) = {result}")

    except ValueError:
        print("Error: Invalid input. Please enter an integer.", file=sys.stderr)
    except Exception as e:
        print(f"An error occurred: {e}", file=sys.stderr)

if __name__ == "__main__":
    solve()
