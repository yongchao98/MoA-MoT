import math

def solve():
    """
    Calculates the number of iterations for the Schwarz Relaxation Method to converge
    for the 1D wave equation.
    """
    # --- Problem Parameters ---
    # c: wave propagation speed
    # T: The time instant up to which convergence is required
    # M: The size of the overlap region (M = b - a)
    
    # You can change these values to match a specific problem
    c = 1.0
    T = 20.0
    M = 0.8

    # The formula for the number of iterations (N) for the parallel Schwarz method is:
    # N = ceil(c * T / M)
    
    # Calculate the argument for the ceiling function
    arg = (c * T) / M
    
    # Calculate the number of iterations
    # The iteration counter starts at 0, so if N iterations are needed, the counter
    # will run from 0 to N-1. The number of iterations is N itself.
    iterations = math.ceil(arg)

    print("For the 1D wave equation, the Schwarz method with absorbing boundary conditions")
    print("converges up to time T after N iterations, where N is determined by the")
    print("wave speed (c), the target time (T), and the overlap size (M).\n")
    print("The governing formula is: N = ceil(c * T / M)\n")
    
    print("Given values:")
    print(f"  Wave speed (c) = {c}")
    print(f"  Target time (T) = {T}")
    print(f"  Overlap size (M) = {M}\n")

    # Print the equation with the numbers plugged in, as requested.
    print("Calculation:")
    print(f"N = ceil({c} * {T} / {M})")
    print(f"N = ceil({arg})")
    print(f"N = {int(iterations)}\n")
    
    print(f"The method needs {int(iterations)} iterations to converge up to time T={T}.")
    print(f"<<<{int(iterations)}>>>")

solve()