import math

# This script demonstrates the principle of option C: using body bias to
# increase the threshold voltage (Vt) and dramatically reduce subthreshold current,
# which in turn increases the pseudo-resistance.

# --- Device & Operating Parameters ---
Vt0 = 0.45  # Nominal threshold voltage in Volts
n = 1.5  # Subthreshold slope factor
Ut = 0.026  # Thermal voltage in Volts (~room temperature)
Vgs = 0.42  # A plausible gate-source voltage for subthreshold operation (must be < Vt)

# --- Body Effect Parameters ---
gamma = 0.4  # Body effect coefficient in V^0.5
phi_f_2 = 0.7  # 2 * Fermi potential in Volts
Vbs_operate = -0.3  # Reverse body-source bias applied during 'operate' phase in Volts

# --- Calculation ---

# 1. Calculate the change in threshold voltage (delta_Vt) due to the body effect.
# The formula is: dVt = gamma * (sqrt(2*phi_f - Vbs) - sqrt(2*phi_f))
delta_Vt = gamma * (math.sqrt(phi_f_2 - Vbs_operate) - math.sqrt(phi_f_2))

# 2. Calculate the new, higher threshold voltage during the operate phase.
Vt_operate = Vt0 + delta_Vt

# 3. The subthreshold current is proportional to: exp((Vgs - Vt) / (n * Ut))
# We can find the ratio of currents (and thus the inverse ratio of resistances)
# without knowing the pre-exponential factor I0.
# Resistance is proportional to 1 / Current.
# Resistance_Ratio = Resistance_with_bias / Resistance_without_bias
#                  = Current_without_bias / Current_with_bias

exponent_no_bias = (Vgs - Vt0) / (n * Ut)
exponent_with_bias = (Vgs - Vt_operate) / (n * Ut)

# The ratio simplifies to exp((Vt_operate - Vt0) / (n * Ut))
resistance_increase_factor = math.exp((Vt_operate - Vt0) / (n * Ut))


# --- Output the results ---
print("--- Analysis for Option C: Dynamic Body Biasing ---")
print(f"Nominal Threshold Voltage (Vt0): {Vt0:.3f} V")
print(f"Applied Reverse Body Bias (Vbs): {Vbs_operate:.3f} V")
print("\nCalculating the increase in threshold voltage (ΔVt):")
print(f"ΔVt = γ * (sqrt(2φ_f - Vbs) - sqrt(2φ_f))")
print(f"ΔVt = {gamma} * (sqrt({phi_f_2} - ({Vbs_operate})) - sqrt({phi_f_2}))")
print(f"ΔVt = {delta_Vt:.3f} V\n")

print("Calculating the new effective threshold voltage (Vt'):")
print(f"Vt' = Vt0 + ΔVt")
print(f"Vt' = {Vt0:.3f} + {delta_Vt:.3f} = {Vt_operate:.3f} V\n")

print("The subthreshold resistance is inversely proportional to the subthreshold current.")
print("The factor by which resistance increases is given by:")
print("Factor = exp(ΔVt / (n * Ut))")
print(f"Factor = exp({delta_Vt:.3f} / ({n} * {Ut}))")
print(f"Factor = exp({delta_Vt:.3f} / {n*Ut:.3f})")
print(f"\nResult: The resistance is increased by a factor of {resistance_increase_factor:.2f}.")
print("\nThis demonstrates that applying a body bias is a highly effective way to increase")
print("the pseudo-resistance, validating strategy C.")
