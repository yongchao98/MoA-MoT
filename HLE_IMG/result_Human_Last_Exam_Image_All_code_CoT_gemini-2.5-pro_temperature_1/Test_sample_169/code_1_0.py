import cmath
import math

# --- Step 1: Define Given Parameters ---
V_target_kV = 0.575  # Target voltage at Bus-W in kV
Z_WF_ohm = 0.01 + 0.05j  # Transmission impedance in Ohms
R_F_ohm = 0.1  # Fault resistance in Ohms
Q_max_MVAR = 10.0  # Maximum reactive power capacity of the compensator in MVAR
PF_min = 0.95  # Minimum power factor (lagging)
harmonic_loss_increase = 0.06  # 6% increase in losses due to harmonics

# --- Step 2: Calculate Equivalent Impedance ---
# The equivalent impedance is the sum of the line impedance and fault resistance
Z_eq_ohm = Z_WF_ohm + R_F_ohm
R_eq_ohm = Z_eq_ohm.real
X_eq_ohm = Z_eq_ohm.imag

# --- Step 3: Determine Power Demand at Bus-W to maintain V_target ---
# We calculate the complex power S = V^2 / Z*
V_target_V = V_target_kV * 1000
V2 = V_target_V**2
Z_eq_conj = Z_eq_ohm.conjugate
S_W = V2 / Z_eq_conj

# Base active power loss (P_loss) and reactive power loss (Q_loss)
P_loss_base_W = S_W.real
Q_loss_VAR = S_W.imag

# --- Step 4: Account for Harmonic Losses ---
# Total active power to be supplied includes harmonic losses
P_total_W = P_loss_base_W * (1 + harmonic_loss_increase)

# --- Step 5: Apply System Constraints to find max Q from wind farm ---
# The power factor constraint applies to the wind farm generator
# PF = P / |S| >= 0.95  => |Q/P| <= tan(acos(PF))
max_tan_phi = math.tan(math.acos(PF_min))
# The wind farm must supply all the active power P_total_W
P_wind_W = P_total_W
# Maximum reactive power the wind farm can supply while meeting the PF constraint
Q_wind_max_VAR = P_wind_W * max_tan_phi

# --- Step 6: Solve for Optimal Reactive Power Injection ---
# The total reactive power demand is Q_loss_VAR.
# Q_loss_VAR = Q_wind + Q_comp
# To minimize Q_comp, we must maximize Q_wind.
# Q_comp_min = Q_loss_VAR - Q_wind_max
Q_opt_VAR = Q_loss_VAR - Q_wind_max_VAR
Q_opt_MVAR = Q_opt_VAR / 1e6

# --- Step 7: Final Verification and Output ---
# Check if the optimal value is within the device's capacity
if 0 <= Q_opt_MVAR <= Q_max_MVAR:
    validity = "This is within the device's capacity."
else:
    validity = "This is outside the device's capacity."

print("--- Nonlinear Optimization Problem Solution ---\n")
print("Objective: Minimize reactive power injection Q_comp.")
print("Constraint: Restore voltage at Bus-W to {:.3f} kV.".format(V_target_kV))
print("Constraint: System power factor >= {:.2f} lagging.".format(PF_min))
print("Constraint: Harmonic losses increase total active power loss by {}%.".format(harmonic_loss_increase * 100))
print("\n--- Calculation Steps ---\n")
print("1. Total equivalent impedance Z_eq = ({:.2f} + j{:.2f}) Ohms".format(R_eq_ohm, X_eq_ohm))

# Converting to MVA for clearer presentation
P_total_MW = P_total_W / 1e6
Q_loss_MVAR = Q_loss_VAR / 1e6
Q_wind_max_MVAR = Q_wind_max_VAR / 1e6

print("2. Total active power required (including harmonic losses): {:.4f} MW".format(P_total_MW))
print("3. Total reactive power required to support voltage: {:.4f} MVAR".format(Q_loss_MVAR))
print("4. Max reactive power from wind farm (at PF={:.2f}): {:.4f} MVAR".format(PF_min, Q_wind_max_MVAR))

print("\n--- Final Equation and Result ---\n")
print("The optimal reactive power Q_opt is found by the equation:")
print("Q_opt = Q_required - Q_wind_max\n")
print("Substituting the calculated values (in MVAR):")
print("{:.4f} = {:.4f} - {:.4f}".format(Q_opt_MVAR, Q_loss_MVAR, Q_wind_max_MVAR))

print("\n----------------------------------------------------")
print("Optimal Reactive Power Injection (Q_opt): {:.4f} MVAR".format(Q_opt_MVAR))
print(validity)
print("----------------------------------------------------")

<<<0.2648>>>