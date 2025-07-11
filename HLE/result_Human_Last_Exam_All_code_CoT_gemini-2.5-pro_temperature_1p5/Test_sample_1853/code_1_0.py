def calculate_gate_capacitance():
    """
    This script derives the formula for the gate capacitance (C_g) based on
    the provided quantum Hall effect measurements.
    """
    # Define symbolic variables for the explanation
    V1 = "V_1"
    B = "B"
    e = "e"
    h = "h"

    # Step 1: Explain the physical model and relationship between V_bg and nu.
    print("Step 1: Relate Gate Voltage and Filling Factor.")
    print("The electron density 'n' is given by n = C_g * V_bg / e and n = nu * e*B/h.")
    print(f"This implies that the gate voltage is proportional to the filling factor: V_bg ∝ nu.\n")

    # Step 2: Identify the specific filling factors from the experimental data.
    print("Step 2: Identify the Filling Factors.")
    print("The observed gate voltages V_1, 3*V_1, and 5*V_1 are in a 1:3:5 ratio.")
    print("Given a total degeneracy of g=4 (spin=2, valley=2), this ratio matches the anomalous Quantum Hall sequence for monolayer graphene (nu = 2, 6, 10, ...).")
    nu_1 = 2
    nu_2 = 6
    nu_3 = 10
    print(f"Thus, the observed levels correspond to filling factors nu = {nu_1}, {nu_2}, and {nu_3}.\n")

    # Step 3: Calculate the change in voltage and filling factor between consecutive levels.
    print("Step 3: Calculate Differences Between Consecutive Observations.")
    # For voltages V_1 and 3*V_1, the corresponding filling factors are nu=2 and nu=6.
    delta_V_coeff = 3 - 1
    delta_nu_val = nu_2 - nu_1
    print(f"The change in voltage: Delta_V = 3*{V1} - 1*{V1} = {delta_V_coeff}*{V1}")
    print(f"The change in filling factor: Delta_nu = {nu_2} - {nu_1} = {delta_nu_val}\n")

    # Step 4: Derive the gate capacitance formula by equating expressions for density change.
    print("Step 4: Derive the Gate Capacitance (C_g).")
    print("Equating the electrostatic and quantum Hall expressions for a change in density:")
    print("C_g * Delta_V / e = Delta_nu * e * B / h")
    print("Solving for C_g gives: C_g = (Delta_nu / Delta_V) * (e^2 * B / h)\n")

    # Step 5: Substitute the numbers and simplify to get the final result.
    print("Step 5: Substitute Values and Calculate Final Expression.")
    print(f"C_g = ({delta_nu_val} / ({delta_V_coeff} * {V1})) * ({e}^2 * {B} / {h})")
    final_coeff = delta_nu_val / delta_V_coeff
    print(f"C_g = ({delta_nu_val}/{delta_V_coeff}) * ({e}^2 * {B}) / ({h} * {V1})")

    print("\n--- Final Equation ---")
    final_equation = f"C_g = {int(final_coeff)} * {e}^2 * {B} / ({h} * {V1})"
    print(final_equation)

if __name__ == '__main__':
    calculate_gate_capacitance()