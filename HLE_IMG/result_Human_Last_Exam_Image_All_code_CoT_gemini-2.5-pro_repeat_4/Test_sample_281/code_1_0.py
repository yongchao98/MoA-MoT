def solve_conductance():
    """
    This function calculates and displays the formula for the four-terminal
    conductance G_12,34 of the specified quantum Hall device.
    """
    # The variables M and N are treated as symbolic characters.
    M = 'M'
    N = 'N'

    print("Based on the Landauer-Büttiker formalism for the given quantum Hall setup, the calculation is as follows:")
    
    # Derivation summary
    print("\n1. Let the current I be sourced from terminal 1 to 2. The voltage is measured between 3 and 4.")
    print("2. Set terminal 2 as ground, so V_2 = 0. Since terminal 4 is a floating voltage probe connected by an edge channel from terminal 2, V_4 = 0.")
    print("3. The QPC scatters M channels coming from the left (at potential V_1) and M channels from the right (at potential V_4 = 0).")
    print("4. N channels are reflected, and M-N channels are transmitted.")
    print("5. The potential at terminal 3, V_3, is determined by the N reflected channels from the left (at V_1) and M-N transmitted channels from the right (at V_4=0).")
    print(f"   V_3 = (N * V_1 + (M - N) * V_4) / M = ({N} / {M}) * V_1")
    print("6. The measured voltage is therefore V_34 = V_3 - V_4 = (N / M) * V_1.")
    print("7. The source current I can be expressed as I = (e^2/h) * (M-N) * V_1.")
    print("8. The four-terminal conductance G_12,34 is I / V_34.")

    # Final formula
    print("\nThe final formula for the conductance is:")
    print(f"G_12,34 = [ (e^2/h) * ({M} - {N}) * V_1 ] / [ ({N} / {M}) * V_1 ]")
    print("After simplifying by canceling V_1, we get:")
    print(f"G_12,34 = (e^2/h) * ({M} * ({M} - {N})) / {N}")

    print("\nIn this final equation:")
    print(f"  {M}: The total number of spin-degenerate edge states.")
    print(f"  {N}: The number of edge states reflected by the QPC.")
    print(f"  e^2/h: The quantum of conductance.")

if __name__ == "__main__":
    solve_conductance()