import math

def solve_schwarz_iterations():
    """
    Calculates the number of iterations for the Schwarz method to converge
    for the 1D wave equation up to a given time T.
    """
    # Define the parameters of the problem.
    # You can change these values to solve for a different case.
    c = 100.0  # Propagation speed
    T = 0.5    # Final time
    a = 40.0   # Left boundary of the right subdomain
    b = 45.0   # Right boundary of the left subdomain

    # The iteration counter starts from 0. The number of iterations is the
    # final value of the counter + 1. The formula calculates this total number.

    # Calculate the size of the overlap, M
    M = b - a

    # Check for invalid overlap size
    if M <= 0:
        print("Error: The overlap size M (b - a) must be positive.")
        return

    # Calculate the expression c * T / M
    expression_value = (c * T) / M

    # The number of iterations N is the ceiling of the expression
    num_iterations = math.ceil(expression_value)

    # Print the explanation and the detailed calculation
    print("The number of iterations (N) for the Schwarz method to converge is given by the formula:")
    print("N = ceil(c * T / M), where M is the overlap size (b - a).\n")

    print("Calculation Steps:")
    print(f"1. Given parameters:")
    print(f"   Propagation speed (c) = {c}")
    print(f"   Final time (T) = {T}")
    print(f"   Overlap boundaries: a = {a}, b = {b}\n")

    print(f"2. Calculate the overlap size (M):")
    print(f"   M = b - a = {b} - {a} = {M}\n")

    print(f"3. Calculate the number of iterations (N):")
    print(f"   N = ceil(c * T / M)")
    print(f"   N = ceil({c} * {T} / {M})")
    print(f"   N = ceil({c * T} / {M})")
    print(f"   N = ceil({expression_value})")
    print(f"   N = {int(num_iterations)}")

if __name__ == '__main__':
    solve_schwarz_iterations()