import numpy as np
from scipy.integrate import solve_ivp

def blowup_model_ode(t, u, k):
    """
    A simplified shell model to demonstrate the blow-up mechanism.
    
    The model is d(u_n)/dt = k_n * u_{n-1}^2 + (1+t) * k_n^2 * u_n, which captures:
    1. Energy cascade to shell n from n-1 (k_n * u_{n-1}^2).
    2. Exponential amplification by the anti-diffusion term ((1+t) * k_n^2 * u_n).
    
    Args:
        t (float): Current time.
        u (np.ndarray): Array of amplitudes [u_0, u_1, ..., u_{N-1}].
        k (np.ndarray): Array of wavenumbers [k_0, k_1, ..., k_{N-1}].
        
    Returns:
        np.ndarray: The time derivatives du/dt.
    """
    N = len(u)
    dudt = np.zeros(N)
    
    # We assume the largest scale (n=0) is constant, representing the initial energy reservoir.
    dudt[0] = 0
    
    # For higher shells (n > 0), both effects are active.
    for n in range(1, N):
        # Nonlinear energy transfer from the next lower shell
        transfer_term = k[n] * u[n-1]**2
        # Linear amplification from the anti-diffusion term
        amplification_term = (1 + t) * (k[n]**2) * u[n]
        dudt[n] = transfer_term + amplification_term
        
    return dudt

# --- Simulation Parameters ---
# Number of shells to model
num_shells = 5
# Wavenumbers for each shell, increasing exponentially
wavenumbers = np.array([2.0**n for n in range(num_shells)])
# Initial condition: energy is concentrated in the largest scale (shell 0).
# Other shells have a tiny seed value to start the cascade.
u_initial = np.zeros(num_shells)
u_initial[0] = 0.1
u_initial[1:] = 1e-12

# Time interval for the simulation. We choose a short interval
# because blow-up is expected to happen quickly.
t_end = 0.4
t_span = [0, t_end]
t_eval = np.linspace(t_span[0], t_span[1], 101)

# --- Solve the ODE system ---
# We use a robust solver from SciPy. If the solution grows too fast (blows up),
# the solver will struggle and stop before reaching t_end.
solution = solve_ivp(
    fun=blowup_model_ode,
    t_span=t_span,
    y0=u_initial,
    args=(wavenumbers,),
    dense_output=True,
    t_eval=t_eval,
    method='RK45'
)

# --- Print and Analyze Results ---
final_time = solution.t[-1]
final_amplitudes = solution.y[:, -1]

print("--- Shell Model Simulation of the Equation ---")
print(f"Simulation ran from t=0 to t={final_time:.4f}")

if final_time < t_end * 0.99:
    print("\n>>> WARNING: The solver stopped before the specified end time.")
    print("This is a strong numerical indication of a finite-time blow-up, as the solution was growing too rapidly for the solver to handle.")
else:
    print("\nSimulation completed for the full time interval.")

print("\nAmplitudes 'u_n' at the final time:")
for i in range(num_shells):
    print(f"  u_{i} (for k={wavenumbers[i]:.1f}): {final_amplitudes[i]:.3e}")

# Calculate and print the total energy E = sum(u_n^2) at different times
total_energy = np.sum(solution.y**2, axis=0)
print("\nGrowth of Total Energy (Σ u_n^2):")
# Print energy at 5 points in time during the simulation
indices_to_print = np.linspace(0, len(solution.t) - 1, 5, dtype=int)
for i in indices_to_print:
    t_i = solution.t[i]
    E_i = total_energy[i]
    print(f"  t = {t_i:.4f}, Total Energy = {E_i:.3e}")
