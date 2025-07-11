import numpy as np
from scipy.integrate import solve_ivp

def main():
    """
    This script estimates the area of a set of initial conditions that lead to
    a specific blow-up behavior in a system of ODEs using a Monte Carlo method.
    """

    # 1. Define the system of ODEs
    def ode_system(t, y):
        """
        Represents the system of ordinary differential equations.
        y[0] = a(t), y[1] = b(t)
        """
        a, b = y
        da_dt = -b * a
        db_dt = -0.5 * b**2 - a**2 + 6 * (a - 1)
        return [da_dt, db_dt]

    # 2. Set up Monte Carlo simulation parameters
    num_samples = 20000
    a_min, a_max = -1.0, 1.0
    b_min, b_max = 2.0, 3.0
    
    # Integration time span
    t_span = [0, 20.0]
    
    # 3. Define the event function for detecting blow-up
    blowup_threshold = 1000.0
    def blowup_event(t, y):
        # Event triggers when the absolute value of a or b exceeds the threshold.
        # This is a proxy for the solution escaping to infinity.
        return max(abs(y[0]), abs(y[1])) - blowup_threshold
    blowup_event.terminal = True  # Stop integration when the event occurs

    # 4. Run the simulation
    blowup_count = 0
    
    # Seed for reproducibility
    np.random.seed(0)

    for i in range(num_samples):
        # Generate a random initial condition from the entire domain R
        a0 = np.random.uniform(a_min, a_max)
        b0 = np.random.uniform(b_min, b_max)
        y0 = [a0, b0]

        # Solve the ODE system
        sol = solve_ivp(
            fun=ode_system,
            t_span=t_span,
            y0=y0,
            method='DOP853',
            events=blowup_event,
            dense_output=True # Required to get state at event time
        )

        # Check if the integration was stopped by our blow-up event
        if sol.status == 1 and sol.t_events[0].size > 0:
            # An event was triggered. Get the state at the event time.
            final_state = sol.y_events[0][0, :]
            a_final, b_final = final_state

            # Check if it's the specific blow-up type: a -> +inf, b -> -inf
            # The event is triggered when |a| or |b| hits the threshold.
            # To match the condition, a must be positive and b must be negative.
            if a_final > 0 and b_final < 0:
                blowup_count += 1
    
    # 5. Calculate and display the result
    total_area = (a_max - a_min) * (b_max - b_min)
    estimated_area = (blowup_count / num_samples) * total_area

    print("--- Monte Carlo Simulation Results ---")
    print(f"Total number of samples: {num_samples}")
    print(f"Number of initial points leading to specified blow-up: {blowup_count}")
    print(f"Total area of the sampling domain [-1,1]x[2,3]: {total_area:.1f}")
    print("\nFinal estimated area calculation:")
    print(f"m(Omega) ~= ({blowup_count} / {num_samples}) * {total_area:.1f}")
    print(f"m(Omega) ~= {estimated_area:.4f}")

if __name__ == '__main__':
    main()
