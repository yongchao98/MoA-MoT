import math

def solve_schwarz_iterations():
    """
    Calculates the number of iterations for the Schwarz Relaxation Method
    for the 1D wave equation to converge up to a time T.
    """
    # --- User-definable parameters ---
    # The final time instant up to which convergence is required.
    T = 12.5
    # The propagation speed of the wave.
    c = 2.0
    # The size of the overlap region (M = b - a).
    M = 4.0
    # ------------------------------------

    # The formula for the number of iterations N is: N = ceil(T * c / M)
    # We calculate the argument for the ceiling function first.
    ratio = (T * c) / M
    
    # Use math.ceil to find the smallest integer >= ratio.
    # The result is cast to int for a whole number of iterations.
    num_iterations = int(math.ceil(ratio))

    print("Calculation of Schwarz Method Iterations for the 1D Wave Equation")
    print("-" * 60)
    print(f"Given parameters:")
    print(f"  Final Time (T)        = {T}")
    print(f"  Propagation Speed (c) = {c}")
    print(f"  Overlap Size (M)      = {M}")
    print("-" * 60)
    print("The required number of iterations 'N' is the smallest integer satisfying:")
    print("N >= T * c / M")
    print("\nStep-by-step calculation:")
    print(f"1. Substitute the values into the inequality:")
    print(f"   N >= {T} * {c} / {M}")
    print(f"2. Calculate the right-hand side:")
    print(f"   N >= {ratio}")
    print(f"3. Find the smallest integer N (Ceiling function):")
    print(f"   N = ceil({ratio})")
    print(f"   N = {num_iterations}")
    print("-" * 60)
    print(f"Final Answer: The method needs {num_iterations} iterations to converge up to time T={T}.")

solve_schwarz_iterations()