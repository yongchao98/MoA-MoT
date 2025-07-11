import numpy as np
from scipy.optimize import minimize

# 1. System parameters (in per-unit)
S_base = 100.0  # MVA
V_G = 1.0 + 0j   # Grid voltage (p.u.)
V_W_target = 0.575 # Target voltage magnitude at Bus-W (p.u.)
Z_WF = 0.01 + 0.05j # Impedance from Bus-W to Fault (p.u.)
Z_GF = 0.01 + 0.05j # Impedance from Fault to Grid (p.u.), assuming symmetry
R_F = 0.1 + 0j      # Fault resistance (p.u.)
Q_comp_max = 10.0 / S_base # Max reactive power from compensator (p.u.)
pf_min = 0.95      # Minimum lagging power factor
tan_phi_max = np.tan(np.arccos(pf_min)) # Corresponds to max Q_W / P_W
harmonic_loss_factor = 1.06 # 6% increase in losses due to harmonics

# 2. Calculate Thevenin equivalent impedance
Z_GF_parallel_RF = (Z_GF * R_F) / (Z_GF + R_F)
Z_th = Z_WF + Z_GF_parallel_RF
R_th = Z_th.real
X_th = Z_th.imag

# 3. Define the optimization problem
# The vector of variables 'x' is [P_W, Q_W, Q_comp]
def objective_function(x):
    """Minimize the reactive power from the compensator."""
    Q_comp = x[2]
    return Q_comp

# Define the nonlinear equality constraint: |V_W| must be V_W_target
def voltage_constraint(x):
    """
    This function represents the power flow equation that must be satisfied.
    It's derived from |V_W|^2 = |V_G - I_inj * Z_th|^2.
    It must equal zero for the constraint to be met.
    """
    P_W, Q_W, Q_comp = x
    Q_inj = Q_W + Q_comp # Total reactive power injected at Bus-W

    # Derived from power balance: (|V_W|^2 + P_W*R_th + Q_inj*X_th)^2 + (-P_W*X_th + Q_inj*R_th)^2 = |V_G|^2 * |V_W|^2
    lhs = (V_W_target**2 + P_W * R_th + Q_inj * X_th)**2 + (-P_W * X_th + Q_inj * R_th)**2
    rhs = (abs(V_G) * V_W_target)**2
    return lhs - rhs

# Define the inequality constraint for the power factor
def pf_constraint(x):
    """
    Represents the constraint Q_W <= tan(phi) * P_W.
    For the solver, this is written as tan(phi)*P_W - Q_W >= 0.
    """
    P_W, Q_W, _ = x
    return tan_phi_max * P_W - Q_W

# Set up constraints for the solver
cons = (
    {'type': 'eq', 'fun': voltage_constraint},
    {'type': 'ineq', 'fun': pf_constraint}
)

# Set up bounds for the variables [P_W, Q_W, Q_comp]
# P_W must be >= 0. No upper bound is given, so we allow it to be large.
# Q_W must be >= 0 (lagging PF).
# Q_comp must be between 0 and Q_comp_max.
bounds = (
    (0, None),
    (0, None),
    (0, Q_comp_max)
)

# Initial guess for the variables
initial_guess = [0.5, 0.1, 0.05]

# 4. Solve the optimization problem
solution = minimize(objective_function, initial_guess, method='SLSQP', bounds=bounds, constraints=cons)

# 5. Output results
if solution.success:
    P_W_opt, Q_W_opt, Q_comp_opt = solution.x
    Q_inj_opt = Q_W_opt + Q_comp_opt

    print("--- Optimal Solution Found ---")
    print(f"Optimal Reactive Power Injection (Q_opt): {Q_comp_opt * S_base:.4f} MVAR ({Q_comp_opt:.4f} p.u.)")
    print(f"Wind Farm Active Power (P_W): {P_W_opt * S_base:.4f} MW ({P_W_opt:.4f} p.u.)")
    print(f"Wind Farm Reactive Power (Q_W): {Q_W_opt * S_base:.4f} MVAR ({Q_W_opt:.4f} p.u.)")

    # 6. Calculate system losses for the optimal point
    # Find complex voltages and currents to calculate losses
    # sin(delta_W) = (-P_W*X_th + Q_inj*R_th) / (|V_G|*|V_W|)
    # cos(delta_W) = (|V_W|^2 + P_W*R_th + Q_inj*X_th) / (|V_G|*|V_W|)
    sin_delta = (-P_W_opt * X_th + Q_inj_opt * R_th) / (abs(V_G) * V_W_target)
    cos_delta = (V_W_target**2 + P_W_opt * R_th + Q_inj_opt * X_th) / (abs(V_G) * V_W_target)
    V_W_opt = V_W_target * (cos_delta + 1j * sin_delta)

    # Calculate V_F
    Y_WF = 1 / Z_WF
    Y_GF = 1 / Z_GF
    Y_F = 1 / R_F
    V_F_opt = (V_W_opt * Y_WF + V_G * Y_GF) / (Y_WF + Y_GF + Y_F)

    # Calculate currents
    I_WF = (V_W_opt - V_F_opt) / Z_WF
    I_GF = (V_F_opt - V_G) / Z_GF
    I_F = V_F_opt / R_F

    # Calculate fundamental losses
    P_loss_WF = abs(I_WF)**2 * Z_WF.real
    P_loss_GF = abs(I_GF)**2 * Z_GF.real
    P_loss_F = abs(I_F)**2 * R_F.real
    P_loss_fund = P_loss_WF + P_loss_GF + P_loss_F
    
    # Calculate total losses including harmonic effects
    P_loss_total = P_loss_fund * harmonic_loss_factor
    
    # We must formulate the equation clearly.
    # The final equation we solved is a constrained optimization problem:
    # min Q_comp
    # s.t. (|V_W|^2 + P_W*R_th + (Q_W+Q_comp)*X_th)^2 + (-P_W*X_th + (Q_W+Q_comp)*R_th)^2 - (|V_G|*|V_W|)^2 = 0
    #      Q_W - tan(acos(0.95))*P_W <= 0
    #      0 <= Q_comp <= 0.1
    #      P_W >= 0, Q_W >= 0
    
    # Plugging in numbers:
    # R_th = 0.0347, X_th = 0.0842
    # (|V_W|=0.575, |V_G|=1.0)
    # (0.575^2 + P_W*0.0347 + (Q_W+Q_comp)*0.0842)^2 + (-P_W*0.0842 + (Q_W+Q_comp)*0.0347)^2 - (1.0*0.575)^2 = 0
    
    print("\n--- System State at Optimal Point ---")
    print(f"The optimization problem was solved to satisfy the voltage equation:")
    print(f"({V_W_target:.3f}^2 + P_W*{R_th:.4f} + (Q_W+Q_comp)*{X_th:.4f})^2 + (-P_W*{X_th:.4f} + (Q_W+Q_comp)*{R_th:.4f})^2 = ({abs(V_G):.1f}*{V_W_target:.3f})^2")
    
    print(f"Final Voltage at Bus-W: {abs(V_W_opt):.4f} p.u.")
    print(f"Total Fundamental Losses: {P_loss_fund * S_base:.4f} MW ({P_loss_fund:.4f} p.u.)")
    print(f"Total Losses (including 6% harmonics): {P_loss_total * S_base:.4f} MW ({P_loss_total:.4f} p.u.)")
    
    # Note: The solver finds that the objective can be minimized to 0. 
    # This means the voltage can be restored without the compensator, but it requires
    # a very high, potentially unrealistic active power injection from the wind farm.
    # The calculation shows P_W > P_loss, but the large values highlight the severity of the fault.

else:
    print("--- Optimization Failed ---")
    print(solution.message)
    print("\nCould not find a feasible solution. This may indicate that the target voltage of 0.575 p.u. is too low to be achieved even with maximum reactive power support from the available sources under the given fault conditions.")
