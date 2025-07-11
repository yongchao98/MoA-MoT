import cmath
import math

def solve_reactive_power_optimization():
    """
    Formulates and solves the nonlinear optimization problem for reactive power injection.
    """
    # Step 1: Define parameters based on the problem statement and engineering assumptions
    S_base = 100.0  # MVA

    # Impedances from the problem (assumed p.u.)
    R_WF = 0.01
    X_WF = 0.05
    R_F = 0.1

    # Active Power from Wind Farm (Assumed 80% of base)
    P_W = 0.8  # p.u.

    # Target Voltage at Bus-W (Interpreted as 0.95 p.u.)
    V_W_target = 0.95  # p.u.

    # Harmonic Loss Factor (6% increase)
    loss_factor = 1.06

    # Assumed Grid Thevenin Impedance (for SCC=20*S_base)
    Z_grid = complex(0, 0.05)
    V_grid = 1.0  # p.u.

    # Step 2: Apply harmonic loss factor to resistances
    R_WF_eff = R_WF * loss_factor
    Z_WF_eff = complex(R_WF_eff, X_WF)
    R_F_eff = R_F * loss_factor

    # Step 3: Calculate the Thevenin equivalent of the faulted grid seen from Bus-W
    Z_th_system = Z_grid + Z_WF_eff
    V_th_faulted = V_grid * (R_F_eff / (R_F_eff + Z_th_system))
    Z_th_faulted = (R_F_eff * Z_th_system) / (R_F_eff + Z_th_system)

    R_th_f = Z_th_faulted.real
    X_th_f = Z_th_faulted.imag
    V_th_f_mag = abs(V_th_faulted)
    Z_th_f_mag2 = R_th_f**2 + X_th_f**2

    # Step 4: Formulate the quadratic equation for Q_comp
    # From V^4 + V^2*(-2*P_W*R - V_th^2) - 2*V^2*Q_comp*X + (P_W^2 + Q_comp^2)*Z_th^2 = 0
    V2 = V_W_target**2
    V4 = V_W_target**4

    a = Z_th_f_mag2
    b = -2 * V2 * X_th_f
    c = V4 + V2 * (-2 * P_W * R_th_f - V_th_f_mag**2) + P_W**2 * Z_th_f_mag2

    print("--- Solving Nonlinear Optimization Problem ---")
    print("\nAn equivalent Thevenin model for the faulted grid is calculated.")
    print(f"The relationship between voltage and power injection leads to a quadratic equation for Q_comp:")
    print(f"\nEquation: a*Q_comp^2 + b*Q_comp + c = 0")
    print(f"a = |Z_th'|^2 = {a:.5f}")
    print(f"b = -2*|V_W|^2*X_th' = {b:.5f}")
    print(f"c = |V_W|^4 + |V_W|^2*(-2*P_W*R_th' - |V_th'|^2) + P_W^2*|Z_th'|^2 = {c:.5f}")
    print("\nFinal Equation to Solve:")
    print(f"{a:.5f} * Q_comp^2 + ({b:.5f}) * Q_comp + ({c:.5f}) = 0")

    # Step 5: Solve the quadratic equation
    discriminant = b**2 - 4 * a * c

    if discriminant < 0:
        print("\n--- Result ---")
        print("The discriminant is negative. No real solution for Q_comp exists.")
        print("The target voltage is physically unreachable with reactive power injection alone.")
        final_answer = float('nan')
    else:
        sqrt_discriminant = math.sqrt(discriminant)
        Q_sol1 = (-b + sqrt_discriminant) / (2 * a)
        Q_sol2 = (-b - sqrt_discriminant) / (2 * a)

        # Select the minimum positive solution as the optimal one
        solutions = [s for s in [Q_sol1, Q_sol2] if s >= 0]
        if not solutions:
            Q_opt_pu = float('nan')
            print("\n--- Result ---")
            print("No positive reactive power solution found.")
        else:
            Q_opt_pu = min(solutions)
            Q_opt_mvar = Q_opt_pu * S_base
            print("\n--- Result ---")
            print(f"The optimal reactive power required to restore voltage to {V_W_target} p.u. is:")
            print(f"Q_opt = {Q_opt_mvar:.4f} MVAR")
            final_answer = Q_opt_mvar

# Run the solver
solve_reactive_power_optimization()