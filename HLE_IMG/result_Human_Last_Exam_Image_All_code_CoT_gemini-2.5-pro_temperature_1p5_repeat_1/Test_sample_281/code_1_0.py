import sympy as sp

def calculate_conductance():
    """
    Calculates the four-terminal conductance G_12,34 for a Hall bar with a QPC.
    """
    # Define symbolic variables for the number of channels
    # M: total number of spin-degenerate edge states
    # N: number of reflected edge states at the QPC
    M, N = sp.symbols('M N', integer=True, positive=True)

    # Define symbolic variables for voltages
    V1, V2 = sp.symbols('V1 V2')

    # Announce the fundamental constants, e (elementary charge) and h (Planck's constant)
    # The conductance quantum is g = e^2/h
    g = sp.Symbol('g') # Represents e^2/h

    # --- Calculation Steps ---
    print("Step-by-step calculation of the conductance G_12,34:\n")

    # Step 1: Determine potentials of probes based on edge state propagation
    print("Step 1: Determine probe potentials.")
    # Potential at terminal 4 is determined by the states coming from terminal 2
    V4 = V2
    print("- The edge states from terminal 2 (at potential V2) flow to terminal 4.")
    print("  Since terminal 4 is a voltage probe, it equilibrates to this potential.")
    print(f"  Therefore, V4 = V2\n")

    # Potential at terminal 3 is determined by the mixing at the QPC
    # It receives N channels from terminal 5 (at V1) and M-N channels from terminal 4 (at V2)
    V3 = (N * V1 + (M - N) * V2) / M
    print("- The edge states reaching terminal 3 come from two paths due to the QPC:")
    print(f"  - {N} states reflected from the top edge (originally from terminal 1, at potential V1).")
    print(f"  - {M-N} states transmitted along the bottom edge (originally from terminal 2, at potential V2).")
    print("- Terminal 3 equilibrates to the average potential of these incoming states.")
    print(f"  Therefore, V3 = ({N}*V1 + ({M}-{N})*V2) / {M}\n")


    # Step 2: Calculate the measured voltage difference V_34
    print("Step 2: Calculate the measured voltage difference V_34 = V3 - V4.")
    V34 = sp.simplify(V3 - V4)
    print(f"  V34 = (({N}*V1 + ({M}-{N})*V2) / {M}) - V2")
    print(f"  Simplifying this gives: V34 = {sp.pretty(V34, use_unicode=False)}\n")


    # Step 3: Calculate the total current I
    print("Step 3: Calculate the total current I from terminal 1 to 2.")
    # The net current I is the sum of currents in the channels crossing the device.
    # It's (M-N) channels carrying charge between potentials V1 and V2.
    I = (M - N) * (V1 - V2) # The g=e^2/h factor is included at the end.
    print(f"- The net current flowing between the left and right side is carried by {M-N} channels.")
    print("  Each channel contributes g*(V1-V2) to the current, where g = e^2/h.")
    print(f"  Therefore, I = ({M} - {N}) * g * (V1 - V2)\n")


    # Step 4: Calculate the conductance G_12,34
    print("Step 4: Calculate the final conductance G_12,34 = I / V34.")
    # We define G_12,34 = I/V_34. The g factor is from I.
    G_12_34_expr = sp.simplify(I / V34)
    print("  G_12,34 = [({M} - {N}) * g * (V1 - V2)] / [({N}/{M}) * (V1 - V2)]")
    print("- The (V1 - V2) terms cancel, leading to the final expression.\n")

    # Final Result
    print("--- Final Result ---")
    print("The four-terminal conductance G_12,34 is given by the formula:")
    # Deconstruct the final expression to clearly show each component.
    numerator_M = M
    numerator_M_minus_N = M - N
    denominator_N = N
    
    print(f"G_12,34 = (({numerator_M}) * ({numerator_M_minus_N})) / ({denominator_N}) * (e^2/h)")


if __name__ == '__main__':
    calculate_conductance()