import cmath
import numpy as np

# Step 1: Define system parameters in per-unit (p.u.)
# Assuming the given impedances are already in p.u. on a 100 MVA base.
# Zs = Transmission line impedance between Bus A and Bus B
Zs = 0.02 + 0.10j
# Zf = Fault impedance
Zf = 0.15 + 0.00j

# Voltages in p.u.
V_fault_mag = 0.85  # Given fault voltage at Bus B
V_target_mag = 1.0  # Target restored voltage at Bus B

# Base values
S_base = 100.0  # MVA
V_base = 220.0  # kV

# Other parameters
Q_max_p_u = 50.0 / S_base # 0.5 p.u.
harmonic_loss_factor = 1.04 # 4% increase

# Step 2: Model the system
# Thevenin impedance seen from Bus B towards the source is Zs.
# The parallel impedance of the system and the fault is Zp = Zs || Zf
Zp = (Zs * Zf) / (Zs + Zf)

# Step 3: Estimate fault voltage phase angle for a consistent model.
# We assume a pre-fault voltage of 1.0 p.u. at the Thevenin source to find a plausible angle.
# This calculated magnitude will differ from the given 0.85, but we use the phase.
V_th_assumed = 1.0 + 0j
V_fault_calc = V_th_assumed * Zf / (Zs + Zf)
V_fault_angle = cmath.phase(V_fault_calc)

# Create the complex fault voltage vector using the given magnitude and calculated angle
V_fault = cmath.rect(V_fault_mag, V_fault_angle)

# Step 4: Solve for the optimal reactive power Q_opt
# The voltage at Bus B (Vb) after compensation is Vb = V_fault - j*Q / Vb_conj * Zp
# This implies that for Q to be real, phase(Vb - V_fault) - phase(Zp) + phase(Vb_conj) + 90deg = 0
# A more direct approach is to solve the system for the angle of Vb (delta) that makes Q real.
# The equation for Q is: Q = cmath.sqrt(-1) * (Vb - V_fault) * Vb.conjugate() / Zp
# We require Q to be real, which means phase( j * (Vb - V_fault) * Vb.conjugate() ) = phase(Zp)
# Or: 90 + phase(Vb - V_fault) - delta = phase(Zp)
# We look for a solution for delta (angle of Vb) that satisfies this.

# Let's solve the system numerically for delta that makes Q real.
def find_q(delta):
    Vb = cmath.rect(V_target_mag, delta)
    Q_complex = 1j * (Vb - V_fault) * Vb.conjugate() / Zp
    return Q_complex

# We search for delta where the imaginary part of Q_complex is zero (or minimized).
# This is equivalent to matching the phase as derived above.
# We are looking for the minimum positive real part of Q_complex.
deltas = np.linspace(-np.pi, np.pi, 360)
solutions = []
for d in deltas:
    q_sol = find_q(d)
    # We are looking for a solution where Q is real and positive (injected).
    if abs(q_sol.imag) < 1e-3: # Check if Q is real within a tolerance
        solutions.append(q_sol.real)

# The minimal positive Q is the optimal one
Q_opt_p_u = min([q for q in solutions if q > 0])
# We need the corresponding delta for loss calculations
final_delta = 0
for d in deltas:
    q_sol = find_q(d)
    if abs(q_sol.real - Q_opt_p_u) < 1e-3 and abs(q_sol.imag) < 1e-3:
        final_delta = d
        break
        
Vb_final = cmath.rect(V_target_mag, final_delta)

# Step 5: Calculate Power Losses
# Current injected by STATCOM
I_comp = (1j * Q_opt_p_u / Vb_final).conjugate()
# Current through the fault
I_fault = Vb_final / Zf
# Total line current from the source is the sum of currents leaving the bus
I_line = I_comp + I_fault

# Losses in system impedance Zs (P = |I|^2 * R)
P_loss_s = (abs(I_line)**2) * Zs.real
# Losses in fault impedance Zf
P_loss_f = (abs(I_fault)**2) * Zf.real

# Total fundamental frequency losses
P_loss_fundamental_p_u = P_loss_s + P_loss_f

# Step 6: Account for harmonics
P_loss_total_p_u = P_loss_fundamental_p_u * harmonic_loss_factor

# Step 7: Display results in actual units
Q_opt_MVAR = Q_opt_p_u * S_base
P_loss_total_MW = P_loss_total_p_u * S_base

# Final Output
print("--- Optimization Results ---")
print(f"Minimum required reactive power injection (Q_opt): {Q_opt_p_u:.3f} p.u.")
print(f"This corresponds to Q_opt = {Q_opt_MVAR:.2f} MVAR.")
# Check against STATCOM limit
if Q_opt_MVAR > 50.0:
    print(f"(Note: This required power of {Q_opt_MVAR:.2f} MVAR exceeds the STATCOM's maximum capacity of 50 MVAR, suggesting the goal may not be achievable with the given hardware.)\n")

print("--- System Power Loss Calculation ---")
print(f"Total system real power losses with harmonics: {P_loss_total_p_u:.3f} p.u.")
print(f"This corresponds to P_loss = {P_loss_total_MW:.2f} MW.")
print("\n--- Detailed Calculation Breakdown ---")
print("Equation for reactive power (Q_opt):")
print(f"|V_B|^2 = |V_F + I_comp * Z_p|^2, where V_B is the target voltage, V_F is the fault voltage, and I_comp is the STATCOM current.")
print(f"Substituting I_comp = (j * Q_opt / V_B)* and solving gives:")
print(f"Q_opt = {Q_opt_MVAR:.2f} MVAR")
print("\nEquation for Power Loss (P_loss):")
print("P_loss = (|I_line|^2 * R_s + |I_fault|^2 * R_f) * (1 + Harmonic_Loss_Factor)")
print(f"P_loss = (({abs(I_line):.2f}^2 * {Zs.real:.2f}) + ({abs(I_fault):.2f}^2 * {Zf.real:.2f})) * {harmonic_loss_factor}")
print(f"P_loss = {P_loss_total_MW:.2f} MW")
