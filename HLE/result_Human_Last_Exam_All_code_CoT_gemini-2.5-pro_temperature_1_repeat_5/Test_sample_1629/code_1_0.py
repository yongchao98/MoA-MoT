import numpy as np
from scipy.integrate import solve_ivp

def estimate_omega_measure():
    """
    Estimates the measure of the set Omega by numerically solving the ODE system
    for a grid of initial conditions.
    """
    # Define the system of ODEs
    def ode_system(t, y):
        a, b = y
        dadt = -b * a
        dbdt = -b**2 / 2 - a**2 + 6 * (a - 1)
        return [dadt, dbdt]

    # Define the grid of initial conditions
    # We use an odd number for na to include a=0 in the grid
    na = 101
    nb = 51
    a0_vals = np.linspace(-1, 1, na)
    b0_vals = np.linspace(2, 3, nb)

    # Set up integration parameters
    t_span = [0, 20]  # Time interval for integration
    blowup_threshold = 1e6

    # Define an event to stop integration if the solution grows too large
    def blowup_event(t, y):
        # Event is triggered when max(|a|, |b|) exceeds the threshold
        return max(abs(y[0]), abs(y[1])) - blowup_threshold
    blowup_event.terminal = True  # Stop integration when event occurs
    blowup_event.direction = 1     # Trigger when event function is increasing

    blowup_count = 0
    total_count = na * nb

    # Iterate over the grid of initial conditions
    for a0 in a0_vals:
        for b0 in b0_vals:
            y0 = [a0, b0]
            # Solve the ODE
            sol = solve_ivp(
                ode_system, t_span, y0,
                events=blowup_event, dense_output=True,
                rtol=1e-5, atol=1e-8
            )

            # Check if the solution blew up in the specified manner
            # sol.status == 1 means a terminal event was triggered
            if sol.status == 1:
                # Get the state at the time of the event
                final_state = sol.y_events[0][-1]
                final_a, final_b = final_state[0], final_state[1]
                
                # Check if a -> +inf and b -> -inf
                if final_a > 0 and final_b < 0:
                    blowup_count += 1

    # Calculate the estimated measure of Omega
    total_area = (a0_vals[-1] - a0_vals[0]) * (b0_vals[-1] - b0_vals[0])
    omega_measure = (blowup_count / total_count) * total_area

    print("--- Numerical Estimation of m(Ω) ---")
    print(f"Total grid points: {total_count}")
    print(f"Initial conditions leading to blow-up: {blowup_count}")
    print(f"Total area of initial domain: {total_area}")
    print("\nFinal Equation:")
    print(f"m(Ω) ≈ ({blowup_count} / {total_count}) * {total_area}")
    print(f"m(Ω) ≈ {omega_measure}")

if __name__ == '__main__':
    estimate_omega_measure()