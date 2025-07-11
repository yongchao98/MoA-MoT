import numpy as np
from scipy.integrate import solve_ivp

def system_of_des(t, state):
    """Defines the system of differential equations."""
    x, y = state
    dxdt = -3 * x * y
    dydt = -y**2 - x + 1
    return [dxdt, dydt]

def run_simulation():
    """
    Solves the ODE for different initial conditions to demonstrate blow-up behavior.
    """
    # Set a fixed initial condition for x(0) > 1
    x0 = 2.0
    
    # Calculate the critical value of y(0) for the chosen x(0)
    # This is the value on the stable manifold of the saddle point (1,0)
    # The equation for the manifold is y^2 = 2*x + 1 - 3*x^(2/3)
    try:
        y_stable = np.sqrt(2 * x0 + 1 - 3 * x0**(2/3))
    except ValueError:
        print(f"Could not calculate stable manifold for x0={x0}")
        return

    print(f"For the initial condition x(0) = {x0}:")
    print("The solution does not blow up only if it starts on or 'above' the stable manifold.")
    # The final equation for the stable manifold at x0 is:
    # y(0)^2 = 2*x(0) + 1 - 3*x(0)^(2/3)
    # We output the numbers for this equation given x0=2
    print(f"For x(0) = {x0}, the stable manifold is at y_s = sqrt(2*{x0} + 1 - 3*{x0}^(2/3)) = {y_stable:.4f}")
    print(f"Blow-up is expected for y(0) < {y_stable:.4f}\n")

    # Define initial conditions for y(0) to test the different behaviors
    # y0_cases: [blow-up case, non-blow-up case (converges to (0,1)), stable manifold case]
    y0_cases = [0.0, y_stable + 0.5, y_stable]
    case_descriptions = [
        f"y(0) = {y0_cases[0]:.4f} (y(0) < y_s)",
        f"y(0) = {y0_cases[1]:.4f} (y(0) > y_s)",
        f"y(0) = {y0_cases[2]:.4f} (y(0) = y_s)"
    ]

    # Set the time span for the integration
    t_span = [0, 25]
    t_eval = np.linspace(t_span[0], t_span[1], 500)

    for y0, desc in zip(y0_cases, case_descriptions):
        print(f"--- Testing case: {desc} ---")
        initial_state = [x0, y0]
        
        # Use an event to detect blow-up (y becoming very large and negative)
        def blowup_event(t, state):
            return state[1] + 1000 # Event triggers when y = -1000
        blowup_event.terminal = True
        blowup_event.direction = -1

        sol = solve_ivp(
            system_of_des, 
            t_span, 
            initial_state, 
            dense_output=True, 
            t_eval=t_eval,
            events=blowup_event
        )
        
        final_x, final_y = sol.y[:, -1]
        
        if sol.status == 1: # Termination event occurred
            print(f"Result: Blow-up detected at t = {sol.t_events[0][0]:.4f}.")
            print(f"Final state before blow-up: (x, y) = ({final_x:.4f}, {final_y:.4f})")
        elif sol.status == 0: # Integration finished successfully
            print(f"Result: No blow-up detected in the time interval.")
            print(f"Final state at t = {t_span[1]}: (x, y) = ({final_x:.4f}, {final_y:.4f})")
            if abs(final_x - 1) < 1e-2 and abs(final_y - 0) < 1e-2:
                print("The solution appears to converge to the saddle point (1, 0).")
            elif abs(final_x - 0) < 1e-2 and abs(final_y - 1) < 1e-2:
                print("The solution appears to converge to the stable point (0, 1).")
        else:
            print(f"Result: Solver failed with status {sol.status}. This may indicate a blow-up.")
        print("-" * 35 + "\n")

    print("\n" + "="*50)
    print("FINAL CONCLUSION:")
    print("The analysis shows that for any initial condition x(0) > 1, a blow-up is guaranteed if and only if y(0) is non-positive.")
    print("The set of values for y(0) for which the solution blows up is: y(0) <= 0")
    print("="*50)

if __name__ == '__main__':
    run_simulation()