# Define example values for the variables.
# Ctg: Top gate capacitance per area (e.g., in F/m^2)
# Vtg: Top gate voltage (in V)
# Cbg: Back gate capacitance per area (e.g., in F/m^2)
# Vbg: Back gate voltage (in V)

Ctg = 1.5e-2 # F/m^2
Vtg = 1.0    # V
Cbg = 0.8e-2 # F/m^2
Vbg = 2.5    # V

# The displacement field (D) is given by D = Ctg * Vtg + Cbg * Vbg.
# The unit of D will be in Coulombs per square meter (C/m^2).

D_tg = Ctg * Vtg
D_bg = Cbg * Vbg
D_total = D_tg + D_bg

print("The displacement field (D) is calculated from the total induced charge density.")
print("The formula is: D = Ctg * Vtg + Cbg * Vbg")
print("\nUsing the provided values:")
print(f"D = {Ctg} * {Vtg} + {Cbg} * {Vbg}")
print(f"D = {D_tg} + {D_bg}")
print(f"D = {D_total} C/m^2")

# Final numerical answer
# The line below is for the final answer extraction.
# <<<0.035>>>