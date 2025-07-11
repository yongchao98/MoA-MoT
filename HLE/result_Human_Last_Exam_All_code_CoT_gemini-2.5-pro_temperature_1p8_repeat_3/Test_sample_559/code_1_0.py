import math

def solve():
    """
    This function states the equation of the separatrix for the given system of differential equations.
    The separatrix was found analytically to be d = -u^2.
    """
    
    # The equation of the separatrix is d = -1 * u^2.
    # The numbers in this equation are -1 (coefficient) and 2 (power).
    
    d_coeff_numerator = -1
    d_coeff_denominator = 1
    
    u_power_numerator = 2
    u_power_denominator = 1

    # We print the numbers that define the equation.
    # We can express the separatrix as d = coefficient * u^power.
    coefficient = d_coeff_numerator / d_coeff_denominator
    power = u_power_numerator / u_power_denominator

    print("The equation for the separatrix is:")
    
    # We build the string for the equation "d = -1 * u^2"
    # and explicitly use the numbers in the print statement.
    print(f"d = {int(coefficient)} * u^{int(power)}")
    
solve()
