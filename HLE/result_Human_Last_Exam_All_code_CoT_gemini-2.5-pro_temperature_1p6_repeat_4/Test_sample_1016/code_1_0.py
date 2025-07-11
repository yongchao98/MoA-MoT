import math

def solve_schwarz_iterations():
    """
    Calculates the number of iterations for the Schwarz method to converge
    for the 1D wave equation up to a time T.
    """
    # --- Problem Parameters ---
    # You can change these values to match your specific problem.
    # Speed of wave propagation
    c = 1.0
    # Left boundary of the right subdomain
    a = 4.5
    # Right boundary of the left subdomain
    b = 6.0
    # Target time for convergence
    T = 25.0

    # --- Calculation ---
    
    # The domain is Omega = [0,L]
    # The left subdomain is Omega_1 = [0,b]
    # The right subdomain is Omega_2 = [a,L]
    # The method requires b > a for an overlap to exist.
    
    # Calculate the size of the overlap region
    M = b - a

    if M <= 0:
        print("Error: Invalid domain definition. The left subdomain's right boundary 'b' must be greater than")
        print("the right subdomain's left boundary 'a' for an overlap to exist.")
        return

    # The convergence time increases by 2*M/c per iteration for the serial Schwarz method.
    # We need to find the number of iterations 'k' such that k * (2 * M / c) >= T.
    # This leads to k >= T * c / (2 * M).
    # Since k must be an integer, we take the ceiling of the result.
    
    # Calculate the raw value before applying the ceiling function
    raw_iterations = (T * c) / (2 * M)
    
    # The number of iterations is the ceiling of the raw value.
    num_iterations = math.ceil(raw_iterations)

    # --- Output ---
    
    print("This script calculates the number of iterations 'k' for the Schwarz method to converge.")
    print("The formula is: k = ceil(T * c / (2 * M))")
    print("\n--- Given Parameters ---")
    print(f"Wave propagation speed (c):         {c}")
    print(f"Target convergence time (T):        {T}")
    print(f"Left subdomain boundary (b):        {b}")
    print(f"Right subdomain boundary (a):       {a}")
    
    print("\n--- Calculation Steps ---")
    print(f"1. Calculate the overlap size M = b - a")
    print(f"   M = {b} - {a} = {M}")
    
    print(f"\n2. Use the convergence formula k >= T * c / (2 * M)")
    print(f"   k >= {T} * {c} / (2 * {M})")
    print(f"   k >= {T * c} / {2 * M}")
    print(f"   k >= {raw_iterations}")
    
    print(f"\n3. Since 'k' must be an integer, we take the ceiling of the value.")
    print(f"   k = ceil({raw_iterations})")
    print(f"   k = {int(num_iterations)}")
    
    print("\n-------------------------")
    print(f"Required number of iterations: {int(num_iterations)}")
    print("-------------------------")


solve_schwarz_iterations()