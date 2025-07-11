import math

def calculate_schwarz_iterations():
    """
    Calculates the number of iterations for the Schwarz method to converge
    for the 1D wave equation up to a given time T.
    """
    # Define the parameters for the wave equation and domain decomposition.
    # These are example values.
    c = 343.0  # Propagation speed (e.g., speed of sound in m/s)
    T = 0.1    # Total simulation time to converge up to (s)
    M = 5.0    # Overlap size (m)

    # The number of iterations 'k' is the smallest integer such that k >= (c * T) / M.
    # This is calculated using the ceiling function.

    # 1. Calculate the argument of the ceiling function
    arg = (c * T) / M

    # 2. Apply the ceiling function to get the number of iterations
    num_iterations = math.ceil(arg)

    # Print the explanation and the final result, showing each number.
    print("This script calculates the number of iterations 'k' for the Schwarz method to converge.")
    print("The formula is: k = ceil(c * T / M)\n")
    print(f"Given values:")
    print(f"  Propagation speed (c) = {c} m/s")
    print(f"  Target time (T)       = {T} s")
    print(f"  Overlap size (M)      = {M} m\n")
    print("Calculation steps:")
    # The prompt asks to output each number in the final equation.
    print(f"  k = ceil({c} * {T} / {M})")
    print(f"  k = ceil({c * T} / {M})")
    print(f"  k = ceil({round(arg, 4)})")
    print(f"  k = {int(num_iterations)}\n")
    print(f"The method needs {int(num_iterations)} iterations to converge up to time T = {T}s.")

# Execute the function
calculate_schwarz_iterations()