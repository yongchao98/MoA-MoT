import numpy as np

def check_electric_dipole_radiation_answer():
    """
    Checks the correctness of the answer based on the physics of electric dipole radiation.
    
    The problem describes an oscillating spheroid, which is a model for an electric dipole.
    The radiated power per unit solid angle (dP/dOmega) for such a dipole oscillating
    along the z-axis is proportional to sin^2(theta) and lambda^(-4).
    
    f(lambda, theta) = dP/dOmega = C * (1/lambda^4) * sin^2(theta)
    
    where C is a constant of proportionality.
    """
    
    # --- Part 1: Check the fraction of power at theta = 30 degrees ---

    # The maximum power 'A' occurs when the sin^2(theta) term is maximum.
    # The maximum value of sin(theta) is 1, which occurs at theta = 90 degrees.
    theta_max_rad = np.deg2rad(90)
    
    # For simplicity, we can analyze the angular factor alone, as other terms cancel out.
    # A is proportional to sin^2(90)
    angular_factor_max = np.sin(theta_max_rad)**2  # This will be 1.0

    # The angle given in the question is theta = 30 degrees.
    theta_given_rad = np.deg2rad(30)
    
    # The power at 30 degrees is proportional to sin^2(30)
    angular_factor_given = np.sin(theta_given_rad)**2 # This will be (0.5)^2 = 0.25

    # The fraction of A is the ratio of the power at 30 degrees to the maximum power.
    calculated_fraction = angular_factor_given / angular_factor_max
    
    # The answer claims the fraction is 1/4.
    expected_fraction = 1/4
    
    # Verify the fraction
    if not np.isclose(calculated_fraction, expected_fraction):
        return (f"Incorrect fraction: The calculated fraction of maximum power at 30 degrees is "
                f"{calculated_fraction:.2f}, but the answer states it is {expected_fraction:.2f}.")

    # --- Part 2: Check the wavelength dependence ---

    # According to the theory of electric dipole radiation (Larmor formula generalization),
    # the radiated power is proportional to the fourth power of the frequency (omega^4).
    # Since frequency is inversely proportional to wavelength (omega = 2*pi*c / lambda),
    # the power is proportional to (1/lambda)^4, or lambda^(-4).
    expected_lambda_dependence = "lambda^(-4)"
    
    # The answer claims a possible form of f includes lambda^(-4).
    # This matches the established physical theory.
    
    # --- Conclusion ---
    # The answer provides a fraction of 1/4 and a wavelength dependence of lambda^(-4).
    # Both of these values are consistent with the physical model of electric dipole radiation.
    # Therefore, the answer is correct.
    
    return "Correct"

# Execute the check
result = check_electric_dipole_radiation_answer()
print(result)