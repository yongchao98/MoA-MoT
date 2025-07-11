import numpy as np
from scipy.optimize import minimize

# 1. System Parameters (in per unit)
V_target = 0.575  # Target voltage at Bus-W [p.u.]
V_s = 1.0         # Assumed effective source voltage [p.u.]
R_eq = 0.01       # Equivalent resistance [p.u.]
X_eq = 0.05       # Equivalent reactance [p.u.]
Q_max = 0.1       # Maximum reactive power injection [p.u.]
PF_min = 0.95     # Minimum power factor (lagging)
loss_factor = 1.06 # Harmonic loss factor

# Effective reactance including harmonic losses
X_eq_eff = X_eq * loss_factor

# 2. Objective Function
# We want to minimize Q_W, which is the second variable in our list x = [P_W, Q_W].
def objective(x):
    P_W, Q_W = x
    return Q_W

# 3. Constraint Definitions
def constraints(x):
    P_W, Q_W = x
    
    # Physics Equality Constraint: V_s^2 - 2*(P*R + Q*X) - V_target^2 = 0
    voltage_constraint = V_s**2 - 2 * (P_W * R_eq + Q_W * X_eq_eff) - V_target**2
    
    # Power Factor Inequality Constraint: P_W / |S_W| >= PF_min
    # This is equivalent to Q_W/P_W <= tan(acos(PF_min))
    # Or, PF_min * P_W - Q_W >= 0 if P_W is positive
    # To handle P_W=0, we write it as P_W^2 * (1/PF_min^2 - 1) >= Q_W^2
    pf_tan_sq = np.tan(np.arccos(PF_min))**2
    pf_constraint = pf_tan_sq * P_W**2 - Q_W**2
    
    return voltage_constraint, pf_constraint

# Define the constraints for the solver
# Constraint 1: Voltage equation must be zero
cons_eq = {'type': 'eq', 'fun': lambda x: constraints(x)[0]}
# Constraint 2: Power factor must be >= 0.95 (our formulation makes this >= 0)
cons_ineq = {'type': 'ineq', 'fun': lambda x: constraints(x)[1]}

# Bounds for the variables [P_W, Q_W]
# P_W must be non-negative
# Q_W must be non-negative (lagging PF) and cannot exceed Q_max
bnds = ((0, None), (0, Q_max))

# Initial guess for [P_W, Q_W]
x0 = [0.5, 0.05]

# 4. Solve the optimization problem
solution = minimize(objective, x0, method='SLSQP', bounds=bnds, constraints=[cons_eq, cons_ineq])

# 5. Output the results
if solution.success:
    P_opt, Q_opt = solution.x
    V_final_sq = V_s**2 - 2 * (P_opt * R_eq + Q_opt * X_eq_eff)
    V_final = np.sqrt(V_final_sq)
    PF_final = P_opt / np.sqrt(P_opt**2 + Q_opt**2) if P_opt > 0 else 0
    Q_opt_MVAR = Q_opt * 100 # Convert back to MVAR

    print("--- Optimal Solution ---")
    print(f"Optimal Reactive Power Injection (Q_opt): {Q_opt:.4f} p.u. ({Q_opt_MVAR:.2f} MVAR)")
    print(f"Required Active Power (P_W): {P_opt:.4f} p.u.")
    print("\n--- Verification of Constraints ---")
    print(f"Resulting Bus-W Voltage: {V_final:.4f} p.u. (Target: {V_target})")
    print(f"Resulting Power Factor: {PF_final:.4f} (Constraint: >= {PF_min})")
    print(f"Reactive Power within Limit: {Q_opt <= Q_max} ({Q_opt:.4f} <= {Q_max})")

    print("\n--- Final Equation with Optimal Values ---")
    # V_s^2 - 2 * (P_W * R_eq + Q_W * X_eq_eff) = V_W^2
    print("The final state satisfies the system physics equation:")
    print(f"{V_s:.2f}**2 - 2 * ({P_opt:.4f} * {R_eq:.2f} + {Q_opt:.4f} * {X_eq_eff:.4f}) = {V_final_sq:.4f}")
    print(f"And the target voltage squared is: {V_target**2:.4f}")
    
else:
    print("Optimization failed.")
    print(solution.message)

# Extract final numerical answer for the platform
final_answer = solution.x[1] * 100 if solution.success else -1
# <<<f"{final_answer:.1f} MVAR">>>