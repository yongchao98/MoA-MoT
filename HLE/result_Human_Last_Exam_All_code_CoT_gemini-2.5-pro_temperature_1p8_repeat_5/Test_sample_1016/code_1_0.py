import math

def solve_schwarz_iterations():
    """
    Calculates the number of iterations needed for the Schwarz Relaxation Method
    for the 1D wave equation to converge up to a time T.
    """

    # --- Problem Parameters ---
    # You can modify these values to fit your specific problem.

    # T: The final time up to which the solution must be converged.
    # For example, if you want the solution to be accurate for t from 0 to 12.5s, T = 12.5.
    T = 12.5  # units: seconds

    # c: The propagation speed of the wave in the medium.
    c = 100.0   # units: meters/second

    # M: The size of the overlap region between subdomains (M = b - a).
    M = 50.0   # units: meters

    # --- Calculation ---
    # The convergence theory states that k iterations are sufficient for convergence
    # up to time t = k * (M/c).
    # To be converged up to time T, we need k * (M/c) >= T.
    # This gives the condition: k >= (T * c) / M.
    # Since k must be an integer, we take the ceiling of the expression.

    ratio = (T * c) / M
    iterations = math.ceil(ratio)

    # --- Output ---
    print("To find the required number of iterations (k) for convergence up to time T,")
    print("we use the formula: k = ceil(T * c / M)")
    print("-" * 50)
    print(f"Given the parameters:")
    print(f"  Final Time (T)   = {T} s")
    print(f"  Wave Speed (c)   = {c} m/s")
    print(f"  Overlap Size (M) = {M} m")
    print("-" * 50)
    # The final print statement shows the full equation with numbers, as requested.
    print("The final equation is:")
    print(f"k = ceil({T} * {c} / {M}) = ceil({T * c} / {M}) = ceil({ratio:.4f})")
    print(f"k = {iterations}")
    print("\nConclusion:")
    print(f"The method needs {iterations} iterations to be converged up to time T = {T}.")


# Execute the function to see the result
solve_schwarz_iterations()