import math

def solve_schwarz_iterations():
    """
    Calculates the number of iterations for the Schwarz Relaxation Method
    for the 1D wave equation to converge up to a time T.
    """
    # Define the parameters of the problem
    # L: Total length of the domain
    # a: Start of the right subdomain
    # b: End of the left subdomain
    # c: Propagation speed of the wave
    # T: Final time up to which convergence is required
    L = 20.0
    a = 8.0
    b = 12.0
    c = 2.0
    T = 50.0

    # Calculate the size of the overlap M
    M = b - a
    
    if M <= 0:
        print("Error: The overlap size M must be positive (b > a).")
        return

    # Calculate the number of iterations k using the formula k = ceil(T * c / M)
    # The iteration counter starts at 0, and k is the number of updates performed.
    # After k updates, the solution is correct up to time t = k * (M / c).
    # We need t >= T, so k * M / c >= T.
    iterations = math.ceil(T * c / M)

    # Print the parameters and the result, including the equation
    print("Problem Parameters:")
    print(f"  Domain Length L = {L}")
    print(f"  Subdomain Ω_1 = [0, {b}]")
    print(f"  Subdomain Ω_2 = [{a}, {L}]")
    print(f"  Overlap Size M = b - a = {M}")
    print(f"  Wave Propagation Speed c = {c}")
    print(f"  Target Time T = {T}")
    print("-" * 30)

    # Output the final equation and result
    print("To reach convergence up to time T, the number of iterations k is given by the formula:")
    print("k = ceil(T * c / M)")
    print("\nPlugging in the values, we get:")
    # The final print statement shows each number in the equation as requested
    print(f"k = ceil({T} * {c} / {M}) = {int(iterations)}")


# Execute the function to solve the task
solve_schwarz_iterations()