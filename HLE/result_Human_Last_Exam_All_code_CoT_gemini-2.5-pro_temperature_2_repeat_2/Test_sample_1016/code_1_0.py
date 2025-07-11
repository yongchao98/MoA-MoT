import math

def solve_schwarz_iterations():
    """
    Calculates the number of iterations for the Schwarz Relaxation Method
    for the 1D wave equation to converge up to a time T.
    """
    # --- Problem Parameters ---
    # These values can be changed by the user.

    # Propagation speed of the wave (e.g., speed of sound in m/s)
    c = 343.0

    # Overlap size between subdomains in meters. Must be > 0 for convergence.
    M = 5.0

    # The time up to which the solution needs to be converged, in seconds.
    T = 0.1

    # --- Calculation ---
    print("--- Schwarz Method Convergence Analysis ---")
    print(f"Wave propagation speed (c): {c} m/s")
    print(f"Overlap size (M): {M} m")
    print(f"Desired convergence time (T): {T} s")
    print("-" * 41)

    # A positive overlap M is required for the method to exchange information
    # and converge in finite time for T > 0.
    if M <= 0:
        if T > 0:
            print("\nError: Overlap size M must be positive for convergence.")
        else:
            # If T is also 0, no time has passed, so 0 iterations are needed.
            print("\nSince T = 0, no time evolution is needed.")
            print("Number of iterations required: 0")
        return

    # If T is zero or negative, no iterations are needed.
    if T <= 0:
        print("\nSince T <= 0, no time evolution is needed.")
        print("Number of iterations required: 0")
        return

    # Step 1: Calculate the total number of update rounds (N) needed.
    # N = ceil( (T * c) / (2 * M) )
    numerator = T * c
    denominator = 2 * M
    num_rounds = math.ceil(numerator / denominator)

    # Step 2: Determine the final iteration index k.
    # The counter k runs from 0 to N-1.
    # The final index is max(0, N-1).
    k_final = max(0, num_rounds - 1)

    # --- Output Results ---
    print("\nStep 1: Calculate the number of full update rounds (N) required.")
    print(f"Formula: N = ceil( (T * c) / (2 * M) )")
    print(f"Calculation: N = ceil( ({T} * {c}) / (2 * {M}) )")
    print(f"             N = ceil( {numerator} / {denominator} )")
    print(f"             N = ceil( {numerator / denominator:.4f} )")
    print(f"             N = {num_rounds}")

    print("\nStep 2: Determine the final iteration index (k).")
    print("The iteration counter k starts at 0, so for N rounds, the indices are 0, 1, ..., N-1.")
    print("The final iteration index needed is k = N - 1 (or 0 if N is 0).")
    print(f"Calculation: k = {num_rounds} - 1 = {k_final}")

    print(f"\n<<<The method needs to run up to iteration {k_final} to achieve convergence at time T.>>>")
    
# Execute the function
solve_schwarz_iterations()