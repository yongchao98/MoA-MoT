import numpy as np
from scipy.integrate import solve_ivp
import sys

# Suppress RuntimeWarning from overflow in exp
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)

def solve_task():
    """
    Solves the system of ODEs for a grid of initial conditions
    and estimates the area of the set Omega that leads to blow-up.
    """

    # The system of ordinary differential equations is:
    # b'(t) = -1/2 * b^2(t) - e^t * a^2(t) - a(t)
    # a'(t) = -1 * b(t) * a(t)
    # We output the coefficients as requested.
    print("Coefficients in the ODE system:")
    print(f"Coefficient of b^2(t) in b'(t): {-0.5}")
    print(f"Coefficient of a^2(t) in b'(t): {-1.0} (multiplied by e^t)")
    print(f"Coefficient of a(t) in b'(t): {-1.0}")
    print(f"Coefficient of b(t)a(t) in a'(t): {-1.0}")
    print("-" * 20)


    def ode_system(t, y):
        """Defines the system of ODEs."""
        a, b = y
        # When numbers get very large, np.exp can overflow.
        # This is part of the "blow-up" and we can handle it.
        try:
            exp_t = np.exp(t)
        except OverflowError:
            exp_t = np.inf
        
        dadt = -b * a
        dbdt = -0.5 * b**2 - exp_t * a**2 - a
        return [dadt, dbdt]

    def event_b_zero(t, y):
        """Event function to detect when b(t) crosses zero."""
        return y[1]
    
    # The event is terminal: integration stops when b(t) = 0.
    event_b_zero.terminal = True
    # We are interested in b crossing from positive to negative.
    event_b_zero.direction = -1

    # Number of random points for Monte Carlo estimation
    num_samples = 4000
    
    # We established that blow-up can only occur for a(0) > 0.
    # The domain for these initial conditions is (0, 1] x [10, 20].
    domain_area = (1.0 - 0.0) * (20.0 - 10.0)
    
    # Time span for integration.
    # If b(t) doesn't cross zero by t=10, it's unlikely to do so later.
    t_span = [0, 10]
    
    blow_up_count = 0

    # Set a seed for reproducibility
    np.random.seed(0)
    
    # Generate random initial conditions
    a0_samples = np.random.uniform(1e-6, 1.0, num_samples)
    b0_samples = np.random.uniform(10.0, 20.0, num_samples)

    for i in range(num_samples):
        y0 = [a0_samples[i], b0_samples[i]]
        
        # Solve the ODE
        sol = solve_ivp(
            ode_system, 
            t_span, 
            y0, 
            events=event_b_zero,
            method='RK45'
        )
        
        # Check if a terminal event (b(t)=0) occurred.
        if sol.status == 1 and len(sol.t_events[0]) > 0:
            blow_up_count += 1
            
    # Calculate the area of the blow-up set Omega
    fraction_blow_up = blow_up_count / num_samples
    m_omega_estimate = fraction_blow_up * domain_area
    
    print(f"Total points sampled in (0, 1] x [10, 20]: {num_samples}")
    print(f"Number of initial conditions leading to blow-up: {blow_up_count}")
    print(f"Estimated fraction of the domain leading to blow-up: {fraction_blow_up:.4f}")
    print(f"The area of the domain (0, 1] x [10, 20] is: {domain_area}")
    print(f"Estimated size m(Omega): {m_omega_estimate:.4f}")

solve_task()