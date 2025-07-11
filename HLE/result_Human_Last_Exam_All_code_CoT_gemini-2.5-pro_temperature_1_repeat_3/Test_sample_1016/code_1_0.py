import math

def calculate_schwarz_iterations():
    """
    Calculates the number of iterations for the Schwarz Relaxation Method
    to converge for the 1D wave equation up to a time T.

    The formula used is k = ceil(T * c / (2 * M)), where:
    - k is the number of iterations.
    - T is the final time instant.
    - c is the wave propagation speed.
    - M is the size of the overlap between subdomains.
    """
    # --- User-defined parameters ---
    # You can change these values to match your specific problem.
    
    # Final time instant
    T = 20.0
    
    # Wave propagation speed
    c = 2.0
    
    # Size of the overlap (M = b - a)
    M = 1.5

    # --- Calculation ---
    
    # The method converges after k iterations, where k is the smallest integer
    # such that k * (2 * M / c) >= T.
    # This gives k >= (T * c) / (2 * M).
    
    numerator = T * c
    denominator = 2 * M
    
    # We need the smallest integer k, so we use the ceiling function.
    # The iteration counter starts at 0, and k=1 is the first iteration.
    num_iterations = math.ceil(numerator / denominator)

    # --- Output the results ---
    print("Problem Parameters:")
    print(f"  Final Time (T)        = {T}")
    print(f"  Wave Speed (c)        = {c}")
    print(f"  Overlap Size (M)      = {M}")
    print("-" * 30)
    
    print("The number of iterations k is calculated using the formula: k = ceil(T * c / (2 * M))")
    print("\nPlugging in the values:")
    print(f"k = ceil({T} * {c} / (2 * {M}))")
    print(f"k = ceil({numerator} / {denominator})")
    print(f"k = ceil({numerator / denominator})")
    print("-" * 30)
    print(f"Number of iterations required: {int(num_iterations)}")

# Execute the function
calculate_schwarz_iterations()