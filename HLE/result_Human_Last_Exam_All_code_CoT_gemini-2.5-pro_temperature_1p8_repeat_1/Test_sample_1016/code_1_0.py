import math

def solve_schwarz_iterations():
    """
    Calculates the number of iterations for the Schwarz method to converge
    for the 1D wave equation.
    """
    # --- Parameters of the problem ---
    # T: The final time instant up to which convergence is required.
    T = 10.0
    # c: The propagation speed of the wave.
    c = 1.0
    # M: The size of the overlap between subdomains (M = b - a).
    M = 1.5

    # --- Core Calculation ---
    # The convergence formula is derived from the time it takes for information
    # to propagate across the overlap region in each iteration.
    # To be converged up to time T, we need k iterations where k >= (T * c) / M.
    # Since k must be an integer, k = ceil(T * c / M).

    # Calculate the value of the ratio (T * c) / M
    value = (T * c) / M

    # Calculate the final number of iterations using the ceiling function
    iterations = math.ceil(value)

    # --- Final Output ---
    print("To find the number of iterations 'k' for the Schwarz method to converge up to time T, we use the formula:")
    print("k = ceil(T * c / M)\n")

    print("Given values:")
    print(f"  T = {T} (Final time)")
    print(f"  c = {c} (Propagation speed)")
    print(f"  M = {M} (Overlap size)\n")

    print("Plugging the values into the equation:")
    # Here we print the equation with the actual numbers, as requested
    print(f"k = ceil({T} * {c} / {M})")
    print(f"k = ceil({T*c} / {M})")
    print(f"k = ceil({value:.4f})")
    print(f"k = {iterations}\n")

    print(f"Thus, the method needs {iterations} iterations to converge up to time T = {T}.")

if __name__ == "__main__":
    solve_schwarz_iterations()
