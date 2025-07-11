import numpy as np
import cmath

# --- Givens and Assumptions ---
# System base power
S_base = 100.0  # MVA
# Transmission line impedance (assumed pu)
Z_S_complex = 0.02 + 0.10j
# External grid voltage (assumed 1.0 pu at the infinite bus A)
V_A = 1.0  # pu
# Bus B voltage during fault
V_B_fault = 0.85  # pu
# Bus B target voltage after compensation
V_B_target = 1.0  # pu
# Harmonic loss increase factor
harmonic_loss_factor = 1.04

# Extract resistance and reactance
R_S = Z_S_complex.real
X_S = Z_S_complex.imag

# --- Step 1: Determine the equivalent fault load (P_FL) ---
# We model the fault as a purely resistive load P_FL that causes V_B to drop to 0.85 pu.
# The relationship between power and voltage gives a quadratic equation for P_FL.
# This equation is derived from the power flow identity sin^2(δ) + cos^2(δ) = 1.
# Equation form: a*P_FL^2 + b*P_FL + c = 0
a1 = R_S**2 + X_S**2
b1 = 2 * R_S * V_B_fault**2 - (2 * R_S * (R_S*0 + X_S*0) / V_B_fault**2) # Simplified for Q=0
c1 = (V_B_fault**4 - (V_A**2 * V_B_fault**2))

# In our specific case (Q_L=0), the equation simplifies. Based on detailed derivation:
# a*x^2+b*x+c=0 for x=P_FL
a_p = (R_S**2 + X_S**2)
b_p = (2 * R_S * V_B_fault**2)
c_p = (V_B_fault**4 - V_A**2 * V_B_fault**2)
# Re-deriving the equation leads to: (R_S^2 + X_S^2)P_FL^2 + (2*R_S*V_B_fault^2)P_FL + (V_B_fault^4 - V_A^2*V_B_fault^2) = 0
# A more robust derivation used was: (X_S^2+R_S^2)P_L^2 + (2R_S V_B^2)P_L + (V_B^4 - V_A^2 V_B^2) = 0
a_pfl = X_S**2 + R_S**2
b_pfl = 2*R_S*V_B_fault**2
c_pfl = V_B_fault**4 - V_A**2*V_B_fault**2
# With numbers: a = 0.0104, b=0.0244, c=-0.2005
a_pfl = 0.0104
b_pfl = 2 * 0.02 * 0.85**2 # = 0.0289
c_pfl = 0.85**4 - 1.0**2 * 0.85**2 #= 0.522 - 0.7225 = -0.2005

print("Step 1: Finding the equivalent fault load P_FL.")
print(f"The quadratic equation for P_FL is: {a_pfl:.4f} * P_FL^2 + {b_pfl:.4f} * P_FL + {c_pfl:.4f} = 0")
p_fl_roots = np.roots([a_pfl, b_pfl, c_pfl])
P_FL = max(p_fl_roots) # Fault must be a positive load
print(f"The equivalent fault load is P_FL = {P_FL:.4f} pu.\n")

# --- Step 2: Determine the optimal reactive power Q_opt ---
# Now find Q_comp needed to bring V_B to 1.0 pu.
# This yields a quadratic equation for Q_comp.
# From detailed derivation: a*Q_c^2 + b*Q_c + c = 0
a_q = R_S**2 + X_S**2
b_q = 2 * X_S * V_B_target**2 - 2 * X_S * P_FL * R_S + 2 * R_S * P_FL * X_S #this term is wrong
b_q_derivation = (2*R_S**2*(-P_FL) + 2*X_S**2*(-P_FL)) # Also wrong.
# Let's use the verified simpler form. (P_B R_S + Q_B X_S + V_B^2)^2 + (Q_B R_S - P_B X_S)^2 = (V_A V_B)^2
P_B_final, V_B_final = P_FL, V_B_target
# (P_B_final*R_S - Q_c*X_S + V_B_final^2)^2 + (-Q_c*R_S - P_B_final*X_S)^2 = (V_A*V_B_final)^2
# Let's expand this for Q_c:
term1_c = V_B_final**2 + P_B_final*R_S # Constant part of first term
term1_q = -X_S # Q_c part of first term
term2_c = -P_B_final*X_S # Constant part of second term
term2_q = -R_S # Q_c part of second term
# (term1_c + term1_q*Q_c)^2 + (term2_c + term2_q*Q_c)^2 = (V_A*V_B_final)^2
# Q_c^2*(term1_q^2+term2_q^2) + Q_c*(2*term1_c*term1_q + 2*term2_c*term2_q) + (term1_c^2+term2_c^2 - (V_A*V_B_final)^2) = 0
a_q = term1_q**2 + term2_q**2
b_q = 2*term1_c*term1_q + 2*term2_c*term2_q
c_q = term1_c**2 + term2_c**2 - (V_A*V_B_final)**2

print("Step 2: Calculating the required reactive power Q_opt.")
print(f"The quadratic equation for Q_comp is: {a_q:.4f} * Q_comp^2 + {b_q:.4f} * Q_comp + {c_q:.4f} = 0")
q_comp_roots = np.roots([a_q, b_q, c_q])
Q_opt_pu = min(q for q in q_comp_roots if np.isreal(q) and q > 0)
Q_opt_MVAR = Q_opt_pu * S_base
print(f"The minimum required reactive power is Q_opt = {Q_opt_pu:.4f} pu ({Q_opt_MVAR:.2f} MVAR).\n")

# --- Step 3: Calculate Total System Losses ---
# Loss is calculated in the final compensated state (V_B=1.0 pu, load = P_FL-jQ_opt)
Q_B_final = -Q_opt_pu

# We need the line current I = (V_A - V_B_final) / Z_S
# First find final angle of V_B
sin_delta_B = (Q_B_final * R_S - P_B_final * X_S) / (V_A * V_B_target)
delta_B_rad = np.arcsin(sin_delta_B)
V_B_final_complex = cmath.rect(V_B_target, delta_B_rad)
# Now calculate current
I_line = (V_A - V_B_final_complex) / Z_S_complex
# Calculate base loss: P = |I|^2 * R
P_loss_base_pu = abs(I_line)**2 * R_S

# Finally, account for the 4% harmonic loss increase
P_loss_total_pu = P_loss_base_pu * harmonic_loss_factor
P_loss_total_MW = P_loss_total_pu * S_base

print("Step 3: Calculating the total system real power losses.")
print(f"The final equation for total losses is: P_loss = |({V_A:.2f} - {V_B_final_complex.real:.3f}{V_B_final_complex.imag:.3f}j) / ({Z_S_complex.real:.2f}+{Z_S_complex.imag:.2f}j)|^2 * {R_S:.2f} * {harmonic_loss_factor:.2f}")
print(f"Total system losses are {P_loss_total_MW:.2f} MW.")
print(f"<<<{P_loss_total_MW:.2f}>>>")