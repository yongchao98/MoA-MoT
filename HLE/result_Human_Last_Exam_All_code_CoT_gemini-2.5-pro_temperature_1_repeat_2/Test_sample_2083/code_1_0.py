import sys

def solve_network_width():
    """
    Calculates the minimum hidden-layer width for a shallow neural network
    to compute the squared norm of an N-dimensional input vector.
    """
    try:
        # Prompt the user for the dimension N
        n_str = input("Enter the dimension of the input vector (N): ")
        N = int(n_str)
        if N <= 0:
            print("Error: N must be a positive integer.", file=sys.stderr)
            return
    except ValueError:
        print(f"Error: Invalid input '{n_str}'. Please enter a positive integer.", file=sys.stderr)
        return
    except (EOFError, KeyboardInterrupt):
        print("\nOperation cancelled by user.", file=sys.stderr)
        return


    # The minimum hidden-layer width H is 2N.
    # This is derived from the fact that approximating an even function like x^2
    # with a non-even activation function like GeLU requires pairing neurons,
    # and the rank of the resulting weight matrix must match the input dimension N.
    # The rank of the matrix is at most H/2, leading to H/2 >= N, or H >= 2N.
    # A constructive proof shows that 2N neurons are also sufficient.
    H = 2 * N

    print(f"For an N-dimensional input vector where N = {N}:")
    print("The minimum required hidden-layer width (H) is given by the equation:")
    # The final print statement outputs each number in the equation as requested.
    print(f"H = 2 * {N}")
    print(f"Result: H = {H}")

if __name__ == "__main__":
    solve_network_width()
