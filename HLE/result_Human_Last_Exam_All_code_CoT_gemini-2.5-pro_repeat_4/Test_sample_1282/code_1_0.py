import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

def enstrophy_ode(t, y, K):
    """
    Defines the simplified ODE for the enstrophy y(t).
    dy/dt = K * y^3 / (1+t)^3
    """
    # To avoid issues with t= -1, although we integrate from t=0
    if 1 + t <= 0:
        return 0
    return K * (y**3) / ((1 + t)**3)

def solve_and_plot(K, y0, t_span, t_eval):
    """
    Solves the ODE and plots the solution.
    """
    sol = solve_ivp(
        fun=enstrophy_ode,
        t_span=t_span,
        y0=[y0],
        args=(K,),
        dense_output=True,
        t_eval=t_eval
    )
    
    T = sol.t
    Y = sol.y[0]
    
    # Check for blow-up (solution grows very large)
    blow_up_time = None
    if np.any(Y > 1e6):
        blow_up_idx = np.where(Y > 1e6)[0][0]
        blow_up_time = T[blow_up_idx]
        
    # --- Analytic solution for verification ---
    # The blow-up occurs if K * y0^2 > 1
    blow_up_condition = K * y0**2
    
    # Denominator of the analytic solution y(t)^2 = y0^2 / (1 - K*y0^2*(1 - 1/(1+t)^2))
    denominator = 1 - blow_up_condition * (1 - 1/(1+T)**2)
    
    # Avoid division by zero or negative sqrt
    denominator[denominator <= 0] = np.nan 
    Y_analytic = y0 / np.sqrt(denominator)

    print(f"Case: K = {K}, y(0) = {y0}")
    print(f"The condition for blow-up in this model is K * y(0)^2 > 1.")
    print(f"For this case, K * y(0)^2 = {blow_up_condition:.4f}")
    
    if blow_up_condition > 1:
        # Calculate theoretical blow-up time
        # T_blowup = sqrt(K*y0^2 / (K*y0^2 - 1)) - 1
        T_blowup_analytic = np.sqrt(blow_up_condition / (blow_up_condition - 1)) - 1
        print(f"The simplified ODE model predicts a blow-up.")
        print(f"Theoretical blow-up time: T = {T_blowup_analytic:.4f}")
        if blow_up_time:
             print(f"Numerical solution appears to blow up around t = {blow_up_time:.4f}")
        else:
             print("Numerical solution did not reach blow-up threshold in the integration time.")
    else:
        print("The simplified ODE model does not predict a blow-up.")

    # --- Print final equation with numbers ---
    # We demonstrate the value at the last time step
    t_final = T[-1]
    y_final = Y[-1]
    # The equation is dy/dt = K * y^3 / (1+t)^3
    # At t_final, the rate is:
    rate_final = K * (y_final**3) / ((1 + t_final)**3)
    print(f"\nFinal state in the simulation at t={t_final:.2f}:")
    print(f"y({t_final:.2f}) = {y_final:.4e}")
    print("The final equation for the rate of change is:")
    print(f"dy/dt = {K:.2f} * ({y_final:.4e})^3 / (1 + {t_final:.2f})^3 = {rate_final:.4e}")

    plt.figure(figsize=(10, 6))
    plt.plot(T, Y, label=f'Numerical Solution (y0={y0})', lw=2)
    plt.plot(T, Y_analytic, 'r--', label=f'Analytic Solution (y0={y0})', lw=2)
    plt.title(f'Solution of the Simplified Enstrophy ODE (K={K})')
    plt.xlabel('Time (t)')
    plt.ylabel('Enstrophy y(t)')
    plt.grid(True)
    plt.legend()
    plt.ylim(0, y0 * 5 if blow_up_condition <= 1 else y0 * 20) # Adjust y-axis for better visualization
    plt.show()


# --- Main execution ---
# Set a value for the constant K
K_const = 0.5 

# Case 1: No blow-up (K * y0^2 <= 1)
# y0^2 <= 1/K => y0 <= sqrt(1/K)
y0_noblowup = np.sqrt(1 / K_const)
solve_and_plot(K=K_const, y0=y0_noblowup, t_span=(0, 10), t_eval=np.linspace(0, 10, 500))

print("\n" + "="*50 + "\n")

# Case 2: Blow-up (K * y0^2 > 1)
y0_blowup = 1.1 * np.sqrt(1 / K_const) # Slightly above the threshold
T_blowup_theory = np.sqrt(K_const * y0_blowup**2 / (K_const * y0_blowup**2 - 1)) - 1
t_end = min(10, T_blowup_theory * 0.99) # Integrate up to near the blow-up time
solve_and_plot(K=K_const, y0=y0_blowup, t_span=(0, t_end), t_eval=np.linspace(0, t_end, 500))
