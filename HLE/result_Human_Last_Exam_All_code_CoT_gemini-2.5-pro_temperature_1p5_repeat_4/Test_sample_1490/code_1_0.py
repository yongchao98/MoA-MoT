import math

def calculate_lcm(numbers):
    """Helper function to calculate LCM of a list of numbers."""
    lcm = 1
    for i in numbers:
        lcm = lcm * i // math.gcd(lcm, i)
    return lcm

def solve_rangoli_problem():
    """
    Solves the rangoli problem based on the provided constraints.
    """
    print("Step 1: Determine the constraints for the number of new curves (D).")
    print(" - For the number of parabolic curves (D/5) to be an integer, D must be a multiple of 5.")
    print(" - For the number of elliptical curves (2D/9) to be an integer, D must be a multiple of 9.")
    print(" - Therefore, D must be a multiple of LCM(5, 9) = 45.")
    print("\nStep 2: Apply the four-way symmetry constraint.")
    print(" - The number of each type of curve must be divisible by 4.")
    print(" - Parabolic curves (D/5) must be divisible by 4, so D must be a multiple of 5 * 4 = 20.")
    print(" - Elliptical curves (2D/9) must be divisible by 4, so D must be divisible by (4 * 9) / 2 = 18.")
    
    # Numbers that the total new curves D must be a multiple of.
    divisors = [20, 18, 9] # From D/5=4k, D/9=2k', and D/9 being an integer
    
    # We just need LCM(20, 18) since 9 is a factor of 18.
    # To be fully rigorous, we find the LCM of all constraints.
    # LCM(D must be mult of 20, D must be mult of 18)
    new_curves_min = calculate_lcm([20, 18])

    print(f"\nStep 3: Calculate the minimum number of new curves that satisfies all constraints.")
    print(f" - D must be a multiple of LCM(20, 18), which is {new_curves_min}.")

    # Calculate the number of curves for each type.
    total_new_curves = new_curves_min
    parabolic_curves = total_new_curves // 5
    elliptical_curves = (total_new_curves * 2) // 9
    circular_curves = total_new_curves - parabolic_curves - elliptical_curves
    
    print(f"\nStep 4: The master must draw {total_new_curves} new curves.")
    print("This total is composed of:")
    print(f" - Parabolic Curves: {total_new_curves} / 5 = {parabolic_curves}")
    print(f" - Elliptical Curves: ({total_new_curves} * 2) / 9 = {elliptical_curves}")
    print(f" - Circular Curves: {total_new_curves} - {parabolic_curves} - {elliptical_curves} = {circular_curves}")

    print("\nFinal Answer: The equation representing the total curves the master must draw is:")
    print(f"{parabolic_curves} + {elliptical_curves} + {circular_curves} = {total_new_curves}")

solve_rangoli_problem()
<<<180>>>