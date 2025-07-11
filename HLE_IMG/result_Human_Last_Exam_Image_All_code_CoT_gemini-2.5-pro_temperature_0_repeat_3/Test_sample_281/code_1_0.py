def calculate_hall_conductance():
    """
    Calculates and prints the four-terminal conductance G_12,34 of a
    quantum Hall device with a quantum point contact (QPC).

    The final derived formula for the conductance is:
    G_12,34 = (M * (M - N) / N) * (e^2/h)
    """
    # --- User-defined parameters ---
    # M: Total number of spin-degenerate edge states.
    #    This corresponds to the filling factor.
    M = 4

    # N: Number of edge states reflected by the QPC.
    #    This value must be an integer between 0 and M (inclusive).
    N = 1
    # -----------------------------

    print("This script calculates the four-terminal conductance G_12,34.")
    print("The formula is derived using the Landauer-Büttiker formalism.")
    print(f"Parameters used: M = {M} (total channels), N = {N} (reflected channels)\n")

    # Input validation
    if not isinstance(M, int) or M < 1:
        print("Error: M must be an integer greater than or equal to 1.")
        return
    if not isinstance(N, int) or not (0 <= N <= M):
        print(f"Error: N must be an integer between 0 and M (inclusive).")
        return

    # The derived formula for the conductance is G = (M * (M - N) / N) * G_0
    # where G_0 = e^2/h is the conductance quantum.

    print("The general formula is: G_12,34 = (M * (M - N) / N) * (e^2/h)")
    print("Substituting the given values into the formula:")
    
    # Handle the special case where N=0 (no reflection)
    if N == 0:
        print(f"G_12,34 = ({M} * ({M} - {N}) / {N}) * (e^2/h)")
        print("\nResult: The conductance is infinite because N=0, leading to division by zero.")
        print("Physically, this means the voltage difference V3 - V4 is zero for a non-zero current.")
    else:
        # Calculate the numerical factor for the conductance
        conductance_factor = M * (M - N) / N
        
        # Print the equation with numbers substituted
        print(f"G_12,34 = ({M} * ({M} - {N}) / {N}) * (e^2/h)")
        print(f"G_12,34 = ({M} * {M - N} / {N}) * (e^2/h)")
        print(f"G_12,34 = {M * (M - N)} / {N} * (e^2/h)")
        print(f"\nResult: G_12,34 = {conductance_factor:.4f} * (e^2/h)")

# Run the calculation function
calculate_hall_conductance()