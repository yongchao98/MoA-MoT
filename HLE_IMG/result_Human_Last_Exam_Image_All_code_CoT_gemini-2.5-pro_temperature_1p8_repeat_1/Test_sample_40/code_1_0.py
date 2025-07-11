import cmath
import math

# Step 1: Given parameters and base values
S_base = 100.0  # MVA
V_B_nominal = 220.0  # kV
V_A = 1.0  # p.u. (Source voltage at Bus A)
V_B_fault = 0.85  # p.u. (Voltage at Bus B during fault)
V_B_final = 1.0  # p.u. (Target voltage at Bus B)

# Assuming Z_S is given in per unit on the 100 MVA base
R_S_pu = 0.02
X_S_pu = 0.10
Z_S_pu = complex(R_S_pu, X_S_pu)

# Step 2: Derive the primary equation from the power-voltage relationship
# (V_final^2 - V_fault^2) / 2 = R_S * Pc + X_S * Qc
# This relates the change in voltage to the STATCOM injection (Pc, Qc)
constant = (V_B_final**2 - V_B_fault**2) / 2

# Step 3: Formulate the second equation based on the loss constraint
# Interpret "4% increase in system losses" as a model for STATCOM losses: Pc = 0.04 * Qc
# This gives a relationship between the real and reactive power of the STATCOM.

# We have a system of two linear equations:
# 1) R_S_pu * Pc_pu + X_S_pu * Qc_pu = constant
# 2) Pc_pu = 0.04 * Qc_pu

# Substitute (2) into (1):
# R_S_pu * (0.04 * Qc_pu) + X_S_pu * Qc_pu = constant
# Qc_pu * (R_S_pu * 0.04 + X_S_pu) = constant

# Step 4: Solve for Q_opt (Qc_pu)
factor = R_S_pu * 0.04 + X_S_pu
Qc_pu = constant / factor
Q_opt_MVAR = Qc_pu * S_base

# Step 5: Calculate the initial system losses
# To do this, we need the fault current. Assume the fault causes a purely reactive power load (Q_L) at Bus B.
# V_fault^2 approx V_A^2 - 2 * (R_S*P_L + X_S*Q_L)
# With P_L = 0, V_fault^2 = V_A^2 - 2 * X_S * Q_L_pu
# 0.85^2 = 1.0^2 - 2 * 0.10 * Q_L_pu
# 0.7225 = 1.0 - 0.2 * Q_L_pu
# 0.2 * Q_L_pu = 1.0 - 0.7225 = 0.2775
Q_L_pu = 0.2775 / 0.2

# Initial fault current I_fault = (S_L / V_fault)*, with S_L = j*Q_L
# I_fault_pu = (complex(0, Q_L_pu) / V_B_fault).conjugate()
S_L_pu_complex = complex(0, Q_L_pu)
I_fault_pu = (S_L_pu_complex / V_B_fault).conjugate()

# Initial line loss (base loss)
P_loss_initial_pu = abs(I_fault_pu)**2 * R_S_pu

# Step 6: Calculate final system losses
# The problem states harmonic effects result in a 4% increase in system losses.
# P_loss_final = 1.04 * P_loss_initial
P_loss_final_pu = 1.04 * P_loss_initial_pu
P_loss_final_MW = P_loss_final_pu * S_base

# Print the results
print("--- Optimization Results ---")
print(f"Objective: Restore Bus B voltage from {V_B_fault*100}% to {V_B_final*100}% of nominal.")

# Using equation R_S * Pc + X_S * Qc = (V_final^2 - V_fault^2)/2
print("\n1. Derivation of required compensation:")
print(f"The voltage restoration requires that:")
print(f"{R_S_pu:.2f} * Pc + {X_S_pu:.2f} * Qc = ({V_B_final:.2f}^2 - {V_B_fault:.2f}^2) / 2")
print(f"{R_S_pu:.2f} * Pc + {X_S_pu:.2f} * Qc = ({V_B_final**2:.4f} - {V_B_fault**2:.4f}) / 2")
print(f"{R_S_pu:.2f} * Pc + {X_S_pu:.2f} * Qc = {constant:.5f}")

print("\n2. Assumption for STATCOM losses:")
print("Assuming STATCOM real power losses (Pc) are 4% of its reactive power generation (Qc):")
print("Pc = 0.04 * Qc")

print("\n3. Calculation of Optimal Reactive Power (Q_opt):")
print(f"Substituting Pc into the equation: {R_S_pu:.2f} * (0.04 * Qc) + {X_S_pu:.2f} * Qc = {constant:.5f}")
print(f"Qc * ({R_S_pu * 0.04:.4f} + {X_S_pu:.2f}) = {constant:.5f}")
print(f"Qc * {factor:.5f} = {constant:.5f}")
print(f"Qc_pu = {constant:.5f} / {factor:.5f} = {Qc_pu:.4f} p.u.")
print(f"Optimal reactive power injection Q_opt = {Qc_pu:.4f} * {S_base:.1f} MVA = {Q_opt_MVAR:.2f} MVAR")
# Note about the capacity limit
print("(Note: This required power of %.2f MVAR exceeds the STATCOM's maximum capacity of 50 MVAR.)" % Q_opt_MVAR)


print("\n4. Calculation of System's Real Power Losses:")
print(f"Initial fault reactive load Q_L_pu = {Q_L_pu:.4f} p.u.")
print(f"Initial fault current |I_fault_pu| = {abs(I_fault_pu):.4f} p.u.")
print(f"Initial system losses P_loss_initial = |I_fault|^2 * R_S = {abs(I_fault_pu):.4f}^2 * {R_S_pu:.2f} = {P_loss_initial_pu:.4f} p.u.")
print("The final system losses are 4% higher than the initial losses.")
print(f"Final system losses = 1.04 * {P_loss_initial_pu:.4f} p.u. = {P_loss_final_pu:.4f} p.u.")
print(f"Final system losses = {P_loss_final_pu:.4f} * {S_base:.1f} MVA = {P_loss_final_MW:.2f} MW")

print("\n--- Final Answer ---")
print(f"The optimal reactive power injection required is Q_opt = {Q_opt_MVAR:.2f} MVAR.")
print(f"The system's final real power losses are {P_loss_final_MW:.2f} MW.")
# Final answer in specified format
# The problem asks for two values. I'll combine them in the final output.
final_answer_str = f"Q_opt = {Q_opt_MVAR:.2f} MVAR, P_loss = {P_loss_final_MW:.2f} MW"
print(f'<<<{final_answer_str}>>>')
