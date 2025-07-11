import math

def calculate_lindhard_function_static_limit():
    """
    Calculates the Lindhard polarization function Pi_0(q=0, omega=0)
    for a 3D electron gas at T=0.

    The calculation is performed in Hartree atomic units (hbar=m_e=e=1)
    for a standard density corresponding to the dimensionless parameter r_s = 1.

    The formula is Pi_0(0,0) = -D(E_F) = -k_F / pi^2, where k_F = (9*pi/4)^(1/3) for r_s=1.
    """
    # Define pi
    pi = math.pi
    
    # The final equation is Pi_0 = - ( (9 * pi / 4)^(1/3) ) / (pi^2)
    # We will output the value of each component of this equation.

    # Calculate the numerator, which is k_F for r_s=1
    numerator_val = (9 * pi / 4)**(1/3)
    
    # Calculate the denominator
    denominator_val = pi**2
    
    # Calculate the final result
    result = -numerator_val / denominator_val
    
    print("The final calculation is based on the equation: Pi_0(0,0) = - (numerator) / (denominator)")
    print(f"The numerator is k_F (for r_s=1 in atomic units): {numerator_val}")
    print(f"The denominator is pi^2: {denominator_val}")
    print(f"The final numerical value is {result}")

if __name__ == '__main__':
    calculate_lindhard_function_static_limit()