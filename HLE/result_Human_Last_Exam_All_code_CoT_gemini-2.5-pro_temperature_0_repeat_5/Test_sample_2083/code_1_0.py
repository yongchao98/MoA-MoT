import sys

def solve_for_width():
    """
    Calculates the minimum hidden-layer width for a shallow neural network
    to compute the squared norm of an N-dimensional vector.
    """
    # The user can change this value to the desired dimension N.
    # If a command-line argument is provided, use it as N.
    if len(sys.argv) > 1 and sys.argv[1].isdigit():
        N = int(sys.argv[1])
    else:
        # Default value for N if not provided
        N = 10

    # The minimum hidden-layer width required is 2*N.
    # This is because the problem decomposes into N independent 1-D problems,
    # each requiring 2 neurons to approximate the non-monotonic x^2 function.
    
    # The two numbers in the equation are the constant '2' and the dimension 'N'.
    constant_factor = 2
    width = constant_factor * N

    print(f"The problem is to find the minimum hidden-layer width (H) for a shallow neural network to compute the squared norm of an N-dimensional vector.")
    print(f"The formula for the minimum width is H = 2 * N.")
    print("-" * 30)
    print(f"Given N = {N}")
    print(f"The final equation is: {constant_factor} * {N} = {width}")
    print("-" * 30)
    print(f"The minimum required hidden-layer width is: {width}")

if __name__ == "__main__":
    solve_for_width()
