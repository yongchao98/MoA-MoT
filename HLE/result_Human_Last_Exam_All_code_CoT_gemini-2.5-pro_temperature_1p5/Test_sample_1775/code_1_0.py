# Define the parameters of the field-effect transistor with example values.
# Ctg: Top gate capacitance per unit area (in Farads per square meter, F/m^2)
Ctg = 1.7e-2

# Vtg: Top gate voltage (in Volts, V)
Vtg = 1.5

# Cbg: Back gate capacitance per unit area (in Farads per square meter, F/m^2)
Cbg = 0.6e-2

# Vbg: Back gate voltage (in Volts, V)
Vbg = 2.0

# The displacement field (D) is the sum of the contributions from the top and back gates.
# The formula is D = Ctg * Vtg + Cbg * Vbg.
# The unit of the displacement field is Coulombs per square meter (C/m^2).

D_total = Ctg * Vtg + Cbg * Vbg

# Print the final equation with all the numbers, as requested.
print("The formula for the total displacement field (D) is: D = Ctg * Vtg + Cbg * Vbg")
print(f"Plugging in the values: D = {Ctg} * {Vtg} + {Cbg} * {Vbg}")
print(f"Result: D = {D_total} C/m^2")

# Final answer for the given example values
# The prompt requires a special format for the final answer.
final_answer = D_total
# <<<final_answer>>>