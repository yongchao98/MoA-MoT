import math

def calculate_frequency_correction_terms():
    """
    Calculates the coefficients of the terms in the nonlinear frequency correction expression.

    The nonlinear frequency correction omega_2 is given by:
    omega_2 = (9*gamma^3 - 4.5*gamma^2 - 3*gamma) / (8 * sqrt(3*gamma))

    This can be expanded into three terms with different powers of gamma.
    """
    
    # Coefficients of the polynomial K = K1*gamma^3 + K2*gamma^2 + K3*gamma
    K1 = 9.0
    K2 = -4.5
    K3 = -3.0
    
    # The denominator is 8 * sqrt(3*gamma)
    # The coefficients of the final expression for omega_2 will be Ki / (8 * sqrt(3))
    
    # Coefficient for the first term (gamma^2.5)
    c1 = K1 / (8 * math.sqrt(3))
    
    # Coefficient for the second term (gamma^1.5)
    c2 = K2 / (8 * math.sqrt(3))
    
    # Coefficient for the third term (gamma^0.5)
    c3 = K3 / (8 * math.sqrt(3))

    print("The nonlinear frequency correction omega_2 is a sum of three terms:")
    print(f"omega_2 = C1 * gamma^2.5 + C2 * gamma^1.5 + C3 * gamma^0.5")
    print("\nCalculated coefficients:")
    print(f"C1 = {K1}/(8*sqrt(3)) = {c1}")
    print(f"C2 = {K2}/(8*sqrt(3)) = {c2}")
    print(f"C3 = {K3}/(8*sqrt(3)) = {c3}")

    print("\nThe 3rd term of the nonlinear correction is the coefficient C3.")
    print(f"The value of the 3rd term is: {c3}")

calculate_frequency_correction_terms()