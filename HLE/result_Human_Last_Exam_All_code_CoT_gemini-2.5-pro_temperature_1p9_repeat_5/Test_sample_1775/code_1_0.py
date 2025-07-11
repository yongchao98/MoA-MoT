def calculate_displacement_field(Ctg, Vtg, Cbg, Vbg):
    """
    Calculates the displacement field in a dual-gate FET.

    Args:
        Ctg (float): Top gate capacitance per unit area (e.g., in F/m^2).
        Vtg (float): Top gate voltage (in V).
        Cbg (float): Bottom gate capacitance per unit area (e.g., in F/m^2).
        Vbg (float): Bottom gate voltage (in V).
    """

    # The total displacement field D is the sum of the contributions from the top and bottom gates.
    # D = D_top + D_bottom = Ctg * Vtg + Cbg * Vbg
    # The dielectric constant of the transistor (epsilon_s) is not needed for this calculation.
    
    # Calculate displacement field
    D_total = Ctg * Vtg + Cbg * Vbg

    # The result has units of charge per area (e.g., C/m^2)

    # Print the equation with all the numbers
    print("The formula for the total displacement field (D) is:")
    print("D = Ctg * Vtg + Cbg * Vbg")
    print("\nPlugging in the values:")
    print(f"D = {Ctg} * {Vtg} + {Cbg} * {Vbg} = {D_total}")

# --- Example Values ---
# You can change these values to match your specific problem.
# Using typical values for demonstration (units are for illustration).
# Capacitance per unit area in microfarads per cm^2 (uF/cm^2)
# Convert to F/m^2 by multiplying by 1e-6 (uF to F) / 1e-4 (cm^2 to m^2) = 1e-2
# Let's use values already in F/m^2 for simplicity.
C_top_gate = 1.5e-2  # F/m^2
V_top_gate = 1.0  # V
C_bottom_gate = 0.5e-2  # F/m^2
V_bottom_gate = 2.0  # V

calculate_displacement_field(C_top_gate, V_top_gate, C_bottom_gate, V_bottom_gate)

# The final numerical answer for the example above
# D = 1.5e-2 * 1.0 + 0.5e-2 * 2.0 = 0.015 + 0.01 = 0.025
final_answer = 1.5e-2 * 1.0 + 0.5e-2 * 2.0
# The requested output format is just the final numerical value.
# However, the user prompt is contradictory. I will just output the number.
# "directly return the answer with the format <<<answer content>>>"
# Let's output this specific example's answer in the required format
# The python code already prints the solution for the user
# So the <<<>>> format needs the equation.
# Re-reading: "Don't include multiple code blocks in one response"
# "remember in the final code you still need to output each number in the final equation!"
# Okay, the code does this.
# "directly return the answer with the format <<<answer content>>> at the end of your response"
# The answer is the formula itself.
final_equation = f"D = {C_top_gate} * {V_top_gate} + {C_bottom_gate} * {V_bottom_gate}"
