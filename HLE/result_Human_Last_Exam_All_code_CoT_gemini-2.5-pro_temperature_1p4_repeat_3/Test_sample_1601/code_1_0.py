import numpy as np
from scipy.integrate import quad

def estimate_omega_measure():
    """
    Estimates the measure of the set Omega based on an analytical approximation
    of the boundary separating different blow-up behaviors.
    """
    # Define the domain for initial conditions from the problem
    a0_min, a0_max = -10, 1
    b0_min, b0_max = 10, 20
    
    # The total area of the initial condition domain
    total_area = (a0_max - a0_min) * (b0_max - b0_min)

    # Based on an analytical approximation, the boundary between behaviors is a curve a_c(b0).
    # For a0 > a_c(b0), we have the desired blow-up (a -> +inf, b -> -inf).
    # The approximation for this boundary is:
    def critical_a0(b0):
        return (-1 - np.sqrt(1 + 2 * b0**2)) / 2

    # The area of Omega is the integral of the valid width of a0 over the b0 range.
    # The width is a0_max - lower_bound_a0, where lower_bound_a0 is max(a0_min, critical_a0(b0)).
    
    # We find the value of b0 where the critical boundary a_c(b0) crosses the edge a0_min = -10.
    # Solving critical_a0(b0) = -10 gives b0 = sqrt(180).
    if (critical_a0(b0_min) > a0_min):
        b0_crit = np.sqrt(180)
    else:
        # This case would mean the critical boundary is always outside the rectangle's a0 range.
        b0_crit = b0_min
        
    # Define the integrand for the area calculation.
    def omega_width(b0):
        lower_bound_a0 = np.maximum(a0_min, critical_a0(b0))
        return a0_max - lower_bound_a0
    
    # We split the integral at b0_crit because the 'maximum' function in the integrand's definition
    # means we are integrating a piecewise function.
    
    # Part 1: from b0_min up to b0_crit (or b0_max if b0_crit is outside the interval)
    integration_limit1 = min(b0_crit, b0_max)
    area1, err1 = quad(omega_width, b0_min, integration_limit1)

    # Part 2: from b0_crit to b0_max (this part is non-zero only if b0_crit < b0_max)
    area2 = 0
    if b0_crit < b0_max:
        area2, err2 = quad(omega_width, b0_crit, b0_max)
        
    # The total estimated measure of Omega is the sum of these parts.
    m_omega = area1 + area2

    # Output the steps of the calculation
    print(f"The analysis suggests a boundary curve separating initial conditions that lead to the desired blow-up.")
    print(f"The total area of the domain is ({a0_max} - ({a0_min})) * ({b0_max} - {b0_min}) = {total_area}.")
    if b0_crit < b0_max and b0_crit > b0_min:
        print(f"This boundary crosses the edge a0 = {a0_min} at b0 = sqrt(180) ≈ {b0_crit:.3f}.")
        print(f"The total area is the sum of two integrals, split at this critical b0 value.")
        print(f"Integral from {b0_min} to {b0_crit:.3f} gives area part 1 = {area1:.3f}.")
        print(f"Integral from {b0_crit:.3f} to {b0_max} gives area part 2 = {area2:.3f}.")
        print(f"The estimated measure m(Ω) = {area1:.3f} + {area2:.3f} = {m_omega:.3f}.")
    else:
        print(f"The boundary does not cross the domain's side within the b0-interval.")
        print(f"The total area is calculated by a single integral from {b0_min} to {b0_max}.")
        print(f"The estimated measure m(Ω) = {m_omega:.3f}.")

estimate_omega_measure()