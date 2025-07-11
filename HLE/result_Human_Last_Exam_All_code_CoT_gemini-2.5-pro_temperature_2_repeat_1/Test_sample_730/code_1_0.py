import math

def solve_parrot_problem():
    """
    Solves the parrot calculation problem by finding an estimation for the rock's mass.
    """
    # 1. State the plan
    print("Yes, the parrot can estimate the mass of the rock with an error of at most 10%.")
    print("Here is the plan to find a calculation using integers as small as possible:")
    print("1. The formula for the mass of the spherical rock is: Mass = Density * (4/3) * pi * Radius^3.")
    print("2. The given values are Density = 0.9 kg/cm3 and Radius = 0.5 cm.")
    print("3. To allow the parrot to calculate this, we must approximate the Density and pi using fractions with small integers.")
    print("4. The volume formula includes (4/3), so the largest integer used must be at least 4. Our goal is to find a valid calculation where the maximum integer is 4.")

    # 2. Define exact values and the target range
    r = 0.5
    rho = 0.9
    true_mass = rho * (4/3) * math.pi * (r**3)
    error_margin = 0.10
    min_mass = true_mass * (1 - error_margin)
    max_mass = true_mass * (1 + error_margin)
    
    # 3. Define the chosen approximations
    # We found that density_approx = 1/1 and pi_approx = 3/1 works.
    rho_num, rho_den = 1, 1
    pi_num, pi_den = 3, 1
    v_const_num, v_const_den = 4, 3
    r_num, r_den = 1, 2

    estimated_mass = (rho_num/rho_den) * (v_const_num/v_const_den) * (pi_num/pi_den) * ((r_num/r_den)**3)
    
    # 4. Present the calculation
    print("\nWe can approximate Density as 1/1 and pi as 3/1.")
    print(f"The estimated mass is {estimated_mass:.3f} kg, which has an error of {abs(estimated_mass/true_mass - 1)*100:.1f}%. This is within the 10% limit.")
    
    print("\nThe calculation for the parrot involves the following numbers:")
    print(f"Density: {rho_num} / {rho_den}")
    print(f"Volume Formula Constant: {v_const_num} / {v_const_den}")
    print(f"Approximation for pi: {pi_num} / {pi_den}")
    print(f"Radius: {r_num} / {r_den}")
    
    print("\nThe final equation for the parrot is:")
    # Using unicode for superscript 3
    print(f"Mass = ({rho_num}/{rho_den}) * ({v_const_num}/{v_const_den}) * ({pi_num}/{pi_den}) * ({r_num}/{r_den})\u00B3")

    # 5. Determine the largest integer
    all_integers = [rho_num, rho_den, v_const_num, v_const_den, pi_num, pi_den, r_num, r_den]
    largest_integer = max(all_integers)
    print(f"\nThe largest integer appearing in this calculation is {largest_integer}.")
    
solve_parrot_problem()

# The final answer format as requested.
print("\n<<<Y4>>>")