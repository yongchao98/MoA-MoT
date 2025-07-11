import numpy as np
from scipy.integrate import solve_ivp

def system_of_odes(t, y):
    """
    Defines the system of ordinary differential equations.
    y[0] corresponds to a(t)
    y[1] corresponds to b(t)
    """
    a, b = y
    a_prime = -b * a
    b_prime = -0.5 * b**2 - a**2 + 6 * (a - 1)
    return [a_prime, b_prime]

def estimate_blowup_set_measure():
    """
    Estimates the measure of the set Omega using a Monte Carlo simulation.
    """
    # Simulation parameters
    N_SAMPLES = 20000  # Number of random initial points to test
    A_RANGE = [-1.0, 1.0]
    B_RANGE = [2.0, 3.0]
    T_MAX = 10.0  # Maximum integration time
    BLOWUP_THRESHOLD = 1.0e3  # Threshold to consider a solution as "blown up"

    # Calculate the total area of the initial domain
    domain_area = (A_RANGE[1] - A_RANGE[0]) * (B_RANGE[1] - B_RANGE[0])

    blowup_count = 0

    # Generate random initial conditions
    # We set a seed for reproducibility of the random numbers
    np.random.seed(42)
    initial_a = np.random.uniform(A_RANGE[0], A_RANGE[1], N_SAMPLES)
    initial_b = np.random.uniform(B_RANGE[0], B_RANGE[1], N_SAMPLES)

    for i in range(N_SAMPLES):
        y0 = [initial_a[i], initial_b[i]]

        # Solve the ODE
        sol = solve_ivp(
            system_of_odes,
            [0, T_MAX],
            y0,
            method='RK45', # A standard solver
            dense_output=True,
            # Stop integration if values get too large to save time
            events=lambda t, y: np.max(np.abs(y)) - BLOWUP_THRESHOLD
        )

        # Check the state at the final time step
        final_a = sol.y[0, -1]
        final_b = sol.y[1, -1]

        # Check for the specified blow-up condition
        # a -> +inf and b -> -inf
        if final_a > BLOWUP_THRESHOLD and final_b < -BLOWUP_THRESHOLD:
            blowup_count += 1
        # Another heuristic for blowup is when solver fails before T_MAX
        # and values are already large and moving in the right direction
        elif sol.status == 1 and sol.t[-1] < T_MAX: # Event triggered
            if final_a > 0 and final_b < 0: # Ensure correct direction
                 # Check velocity at the end to confirm divergence
                v_final = system_of_odes(sol.t[-1], sol.y[:, -1])
                if v_final[0] > 0 and v_final[1] < 0:
                     blowup_count += 1


    # Estimate the measure of Omega
    estimated_measure = (blowup_count / N_SAMPLES) * domain_area

    # Output the results
    print("--- Monte Carlo Simulation Results ---")
    print(f"Total initial points sampled: {N_SAMPLES}")
    print(f"Number of points leading to blow-up (a->+inf, b->-inf): {blowup_count}")
    print(f"Total area of the initial domain [-1,1]x[2,3]: {domain_area}")
    print("\n--- Final Calculation ---")
    print(f"m(Omega) ~= (Points in Omega / Total Points) * Total Area")
    print(f"m(Omega) ~= ({blowup_count} / {N_SAMPLES}) * {domain_area:.1f}")
    print(f"Estimated measure m(Omega): {estimated_measure:.4f}")

if __name__ == '__main__':
    estimate_blowup_set_measure()
