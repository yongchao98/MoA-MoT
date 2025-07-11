import sys

def calculate_conductance(M, N):
    """
    Calculates the four-terminal conductance G_12,34 for a quantum Hall device with a QPC.

    Args:
        M (int): The total number of spin-degenerate edge states.
        N (int): The number of edge states reflected by the QPC.
    """
    # Fundamental constants are treated symbolically, G_0 = e^2/h
    G0_str = "e^2/h"

    print(f"Calculating conductance for M = {M} and N = {N}")
    print("-" * 30)

    # Input validation
    if not isinstance(M, int) or not isinstance(N, int) or M <= 0 or N < 0 or N > M:
        print("Error: Invalid input.")
        print("M must be a positive integer (total channels).")
        print("N must be an integer such that 0 <= N <= M (reflected channels).")
        return

    # Handle special cases
    if N == 0:
        # No reflection, QPC is fully transparent.
        # The voltage difference V3 - V4 becomes zero, leading to infinite conductance.
        print("The QPC is fully transparent (N=0).")
        print("The voltage difference V34 is zero.")
        print(f"G_12,34 = I / (V3 - V4) = infinity")
        return

    if N == M:
        # Full reflection, QPC is fully closed.
        # The current path from 1 to 2 is blocked, leading to zero conductance.
        print("The QPC is fully reflecting (N=M).")
        print("The conductance G_12,34 is zero.")
        print(f"G_12,34 = ({M} * ({M} - {N}) / {N}) * {G0_str} = 0.0 * {G0_str}")
        return

    # General case: 0 < N < M
    # Formula: G = (M * (M - N) / N) * G_0
    conductance_coeff = M * (M - N) / N

    print("The formula for the conductance is:")
    print("G_12,34 = (M * (M - N) / N) * e^2/h")
    print("\nSubstituting the values:")
    # Output each number in the final equation
    print(f"G_12,34 = ({M} * ({M} - {N}) / {N}) * {G0_str}")
    print(f"G_12,34 = ({M} * {M - N} / {N}) * {G0_str}")
    print(f"G_12,34 = {M * (M - N)} / {N} * {G0_str}")
    print(f"G_12,34 = {conductance_coeff:.4f} * {G0_str}")


if __name__ == '__main__':
    # --- User-definable parameters ---
    # M: Total number of spin-degenerate edge states.
    # N: Number of reflected edge states by the QPC.
    # You can change these values to see different results.
    M = 3
    N = 1
    # --- End of user-definable parameters ---

    calculate_conductance(M, N)