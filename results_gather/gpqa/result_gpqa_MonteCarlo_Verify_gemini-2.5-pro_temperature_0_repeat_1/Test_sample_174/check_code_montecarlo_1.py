import math

def check_correctness():
    """
    Checks the correctness of the given answer for the physics problem.

    The problem describes an oscillating spheroid charge distribution, which is a model
    for an oscillating electric dipole. The radiated power per unit solid angle (dP/dΩ)
    for such a dipole is proportional to sin²(θ) / λ⁴, where θ is the angle from the
    oscillation axis and λ is the wavelength.
    """
    
    # The provided answer is 'A'. Let's parse the values from option A.
    # A) 1/4, λ^(-4)
    expected_fraction = 1/4
    expected_lambda_exponent = -4

    # --- Verification Part 1: Fraction of Maximum Power ---

    # The radiated power is proportional to sin²(θ).
    # The maximum power, A, occurs when sin²(θ) is maximum.
    # The maximum value of sin²(θ) is 1, which occurs at θ = 90 degrees (π/2 radians).
    # Let the power be P(θ) = K * sin²(θ) for some constant K.
    # A = P(90°) = K * sin²(90°) = K * 1 = K.
    
    # We need to find the power at θ = 30 degrees.
    theta_deg = 30
    theta_rad = math.radians(theta_deg)
    
    # Power at 30 degrees is P(30°) = K * sin²(30°).
    # sin(30°) = 1/2, so sin²(30°) = (1/2)² = 1/4.
    # P(30°) = K * (1/4).
    
    # The fraction of A radiated at 30 degrees is P(30°) / A.
    calculated_fraction = (0.25) / 1.0

    # Check if the calculated fraction matches the expected fraction.
    if not math.isclose(calculated_fraction, expected_fraction):
        return (f"Incorrect: The fraction of maximum power at 30 degrees is wrong. "
                f"Power is proportional to sin²(θ). Maximum power A is at θ=90° (sin²(90°)=1). "
                f"At θ=30°, the factor is sin²(30°) = (0.5)² = 0.25. "
                f"The fraction should be 0.25/1 = 1/4. The answer states {expected_fraction}.")

    # --- Verification Part 2: Wavelength Dependence ---

    # The theory of electric dipole radiation states that the total radiated power
    # is proportional to the fourth power of the angular frequency (ω⁴).
    # Since ω = 2πc/λ, the power is proportional to (1/λ)⁴ or λ⁻⁴.
    # The power per unit solid angle, f(λ, θ), inherits this dependence.
    calculated_lambda_exponent = -4

    # Check if the calculated exponent matches the expected exponent.
    if calculated_lambda_exponent != expected_lambda_exponent:
        return (f"Incorrect: The dependence on wavelength λ is wrong. "
                f"For electric dipole radiation, power is proportional to ω⁴, where ω=2πc/λ. "
                f"Therefore, power is proportional to (1/λ)⁴ = λ⁻⁴. "
                f"The answer suggests a dependence of λ^({expected_lambda_exponent}).")

    # If both checks pass, the answer is correct.
    return "Correct"

# Run the check
result = check_correctness()
print(result)