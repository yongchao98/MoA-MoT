import math

def calculate_lindhard_function_value():
    """
    Calculates the dimensionless numerical value of the Lindhard polarization function
    at zero frequency and zero momentum transfer for a 3D electron gas.
    
    The Lindhard function in this limit is Pi(0,0) = -g(epsilon_F), where g(epsilon_F)
    is the density of states at the Fermi energy. For a 3D electron gas, this
    is given by g(epsilon_F) = (3 * n) / (2 * epsilon_F), where n is the electron
    density and epsilon_F is the Fermi energy.

    To find a universal numerical value, we evaluate this quantity in the natural
    dimensionless units of the system, where the unit of energy is epsilon_F and
    the unit of volume is 1/n. This makes the parameters n and epsilon_F cancel out.
    """

    # The formula is Pi(0,0) = - (numerator_factor * n) / (denominator_factor * epsilon_F)
    numerator_factor = 3
    denominator_factor = 2

    # In natural units, the calculation for the dimensionless value becomes:
    # Dimensionless Value = - numerator_factor / denominator_factor
    dimensionless_value = -numerator_factor / denominator_factor

    print("The Lindhard polarization function in the static, long-wavelength limit is given by Pi(0,0) = -g(epsilon_F).")
    print("For a 3D electron gas, the density of states at the Fermi level g(epsilon_F) = (3 * n) / (2 * epsilon_F).")
    print("\nTherefore, the expression for the Lindhard function is:")
    print(f"Pi(0,0) = - ({numerator_factor} * n) / ({denominator_factor} * epsilon_F)")
    
    print("\nTo obtain a universal numerical value, we express this in the system's natural units.")
    print("This yields a dimensionless value based on the constant factors in the equation.")
    print("\nThe final equation for the numerical value is:")
    print(f"Value = - {numerator_factor} / {denominator_factor}")
    
    print("\nExecuting the calculation:")
    print(f"The numerical value is: {dimensionless_value}")


# Run the calculation
calculate_lindhard_function_value()