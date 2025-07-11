import numpy as np
from scipy.optimize import minimize

# 1. Define System Parameters in Per-Unit (p.u.)
Z_S_complex = 0.02 + 0.10j
R_S = Z_S_complex.real
X_S = Z_S_complex.imag
Z_S_sq = np.abs(Z_S_complex)**2

Z_F = 0.15 # p.u.
Q_max = 0.5 # p.u. (50 MVAR / 100 MVA)
V_A = 1.0 # p.u.
V_B_fault = 0.85 # p.u.
V_B_target = 1.0 # p.u.
S_base = 100.0 # MVA
pf_min = 0.98
loss_increase_factor = 1.04

# STATCOM power factor constraint tan(phi) = P_comp / Q_comp
tan_phi_max = np.tan(np.arccos(pf_min))

# 2. Formulate the Non-Linear Optimization Problem
# Variables: x = [P_total_fault, Q_total_fault, P_comp, Q_comp]

# 3. Define Objective Function
def objective(x):
    P_total_fault, Q_total_fault, P_comp, Q_comp = x
    return Q_comp

# 4. Establish Constraints
def constraints(x):
    P_total_fault, Q_total_fault, P_comp, Q_comp = x
    
    # Initial Condition Constraint (V_B = 0.85)
    v_b1 = V_B_fault
    eq1 = v_b1**4 + (P_total_fault**2 + Q_total_fault**2)*Z_S_sq + 2*v_b1**2*(P_total_fault*R_S + Q_total_fault*X_S) - (V_A**2 * v_b1**2)

    # Final Condition Constraint (V_B = 1.0)
    # Using a constant impedance load model
    k = (V_B_target / V_B_fault)**2
    P_final = k * P_total_fault - P_comp
    Q_final = k * Q_total_fault - Q_comp
    v_b2 = V_B_target
    eq2 = v_b2**4 + (P_final**2 + Q_final**2)*Z_S_sq + 2*v_b2**2*(P_final*R_S + Q_final*X_S) - (V_A**2 * v_b2**2)
    
    # STATCOM Power Factor Constraint
    # P_comp - tan_phi_max * Q_comp <= 0
    cons_pf = P_comp - tan_phi_max * Q_comp
    
    return [eq1, eq2, cons_pf]

# Bounds for the variables [P_total_fault, Q_total_fault, P_comp, Q_comp]
bnds = ((-20, 20), (-20, 20), (0, None), (0, Q_max))
# Initial guess for the optimizer
x0 = [5.0, 10.0, 0.05, 0.45] 
# Define constraints for the solver
cons = ({'type': 'eq', 'fun': lambda x: constraints(x)[0]},
        {'type': 'eq', 'fun': lambda x: constraints(x)[1]},
        {'type': 'ineq', 'fun': lambda x: -constraints(x)[2]}) # Constraint must be >= 0

# 5. Solve Numerically
solution = minimize(objective, x0, method='SLSQP', bounds=bnds, constraints=cons, options={'disp': False, 'ftol': 1e-9})

if solution.success:
    P_total_fault_opt, Q_total_fault_opt, P_comp_opt, Q_comp_opt = solution.x
    Q_opt_MVAR = Q_comp_opt * S_base

    # 6. Calculate System Losses
    # Determine the final operating point and voltage angle
    k = (V_B_target / V_B_fault)**2
    P_final = k * P_total_fault_opt - P_comp_opt
    Q_final = k * Q_total_fault_opt - Q_comp_opt
    
    # Find cos(delta) and sin(delta) from power flow equations
    M_inv = (1/Z_S_sq) * np.array([[ -R_S, -X_S], [X_S, -R_S]])
    vec_in = np.array([P_final * Z_S_sq / V_B_target, Q_final * Z_S_sq / V_B_target])
    vec_out = M_inv @ vec_in
    A = (vec_out[0] / V_B_target) + V_B_target
    B = vec_out[1] / V_B_target
    cos_delta = A / V_A
    sin_delta = B / V_A
    delta_final = np.arctan2(sin_delta, cos_delta)

    # Calculate final current and losses
    V_A_phasor = V_A
    V_B_phasor = V_B_target * (cos_delta + 1j * sin_delta)
    I_final_phasor = (V_A_phasor - V_B_phasor) / Z_S_complex
    P_loss_base = np.abs(I_final_phasor)**2 * R_S
    P_loss_total = loss_increase_factor * P_loss_base
    P_loss_total_MW = P_loss_total * S_base
    
    # 7. Output the results
    print(f"--- Optimization Results ---")
    print(f"The minimum reactive power injection required from the STATCOM is Q_opt.")
    print(f"Equation: Q_opt = {Q_comp_opt:.4f} * {S_base:.1f} = {Q_opt_MVAR:.2f}")
    print(f"Result: Q_opt = {Q_opt_MVAR:.2f} MVAR")
    
    print(f"\nThe system's real power losses after compensation are calculated.")
    print(f"Equation: P_loss = {loss_increase_factor:.2f} * |I_final|^2 * R_S * S_base")
    print(f"Equation: P_loss = {loss_increase_factor:.2f} * {np.abs(I_final_phasor):.2f}^2 * {R_S:.2f} * {S_base:.1f} = {P_loss_total_MW:.2f}")
    print(f"Result: P_loss = {P_loss_total_MW:.2f} MW")

else:
    print("Optimization failed.")
    print(solution.message)
