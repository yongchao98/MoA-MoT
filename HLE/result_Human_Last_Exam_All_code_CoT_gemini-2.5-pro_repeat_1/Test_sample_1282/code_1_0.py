import numpy as np
from scipy.integrate import solve_ivp

def solve_and_explain():
    """
    Analyzes the blow-up possibility by solving a derived ODE numerically.
    """
    # In the analysis, the problem is reduced to studying the ODE:
    # dY/dτ = K * Y^3 / (1 + 2τ)^2
    # We set the constant K=1.0 for this demonstration.
    K = 1.0

    # Define the ODE system for Y(τ), an upper bound on the squared L2-norm of the vorticity.
    def ode_system(t, y, K_val):
        return [K_val * (y[0]**3) / (1 + 2*t)**2]

    # Event to detect blow-up (when the solution becomes very large)
    def blow_up_event(t, y, K_val):
        return 1e6 - y[0]
    blow_up_event.terminal = True
    blow_up_event.direction = -1

    print("Analyzing the possibility of a finite-time blow-up for the given PDE.")
    print("The analysis leads to a differential inequality for an upper bound Y(τ) on the solution's norm:")
    print("\ndY/dτ ≤ K * Y^3 / (1 + 2τ)^2\n")
    print("We solve the corresponding ODE to demonstrate its behavior for K=1.")
    print("-" * 70)

    # --- Case 1: Small initial data ---
    # The analytical solution of the ODE blows up only if Y(0)^2 > 1/K.
    # So, if Y(0) <= 1/sqrt(K), the estimate remains bounded.
    y0_small_val = 0.95 / np.sqrt(K)
    sol_small = solve_ivp(ode_system, [0, 50], [y0_small_val], args=(K,), method='RK45')
    
    print("Case 1: Small initial data")
    print(f"The analysis predicts the solution remains bounded if Y(0) ≤ {1/np.sqrt(K):.2f}.")
    print(f"Choosing Y(0) = {y0_small_val:.2f}, the numerical solution confirms this.")
    print(f"At τ = {sol_small.t[-1]:.1f}, the solution is Y = {sol_small.y[0, -1]:.4f}, which is finite.")
    print("-" * 70)

    # --- Case 2: Large initial data ---
    # We choose Y(0) such that Y(0)^2 > 1/K, where blow-up of the estimate is possible.
    y0_large_val = 1.5 / np.sqrt(K)
    
    # Calculate the theoretical blow-up time from the analytical solution of the ODE.
    # The formula is T_blowup = 1 / (2 * (K * Y(0)^2 - 1)).
    T_blowup_numerator = 1.0
    T_blowup_denominator = 2.0 * (K * y0_large_val**2 - 1.0)
    T_blowup = T_blowup_numerator / T_blowup_denominator
    
    sol_large = solve_ivp(ode_system, [0, T_blowup * 2], [y0_large_val], args=(K,), method='RK45', events=blow_up_event)

    print("Case 2: Large initial data")
    print(f"The analysis suggests a blow-up is possible if Y(0) > {1/np.sqrt(K):.2f}.")
    print(f"Choosing Y(0) = {y0_large_val:.2f}, the ODE solution should blow up.")
    print("\nThe theoretical blow-up time T_blowup is given by the equation:")
    print(f"T_blowup = {T_blowup_numerator} / (2 * ({K} * {y0_large_val:.2f}^2 - 1.0))")
    print(f"Calculated T_blowup = {T_blowup:.4f}")
    
    if sol_large.status == 1: # Event triggered
        print(f"The numerical solver stopped at τ = {sol_large.t_events[0][0]:.4f}, confirming the blow-up.")
    else:
        print("Numerical solver did not detect a blow-up within the integration interval.")
    print("-" * 70)
    
    print("\nConclusion: Since the standard energy method fails to guarantee global regularity for large data")
    print("and instead predicts a finite-time blow-up for the solution's norm bound, it is considered")
    print("possible for the solution of the PDE to blow up in finite time.")

if __name__ == '__main__':
    solve_and_explain()