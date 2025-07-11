import math

def solve_egg_drop_problem():
    """
    Calculates the time for an egg to fall from a skyscraper using approximations
    with small integers and fractions, suitable for a young learner.
    """
    # Given values from the problem
    distance_from_base = 240  # in meters

    # --- Approximations using small integers/fractions ---
    # The angle is 22.5 degrees. We approximate tan(22.5) ≈ 5/12.
    tan_alpha_num = 5
    tan_alpha_den = 12

    # We approximate gravity g ≈ 10 m/s^2.
    g_approx = 10

    # We will need to approximate sqrt(5). We use sqrt(5) ≈ 9/4.
    sqrt5_approx_num = 9
    sqrt5_approx_den = 4

    # --- Step 1: Calculate the height (h) ---
    print("--- Part 1: Finding the Skyscraper's Height ---")
    print(f"The height 'h' is the distance from the base multiplied by tan(alpha).")
    print(f"We are given the distance = {distance_from_base} m.")
    print(f"We approximate tan(alpha) with the fraction: {tan_alpha_num} / {tan_alpha_den}")
    print(f"Calculation: h = {distance_from_base} * ({tan_alpha_num} / {tan_alpha_den})")
    
    # To keep numbers small, we divide before multiplying
    h_intermediate = distance_from_base // tan_alpha_den
    h = h_intermediate * tan_alpha_num
    print(f"h = ({distance_from_base} / {tan_alpha_den}) * {tan_alpha_num} = {h_intermediate} * {tan_alpha_num} = {h} m")
    print("\n")

    # --- Step 2: Calculate the fall time (t) ---
    print("--- Part 2: Finding the Fall Time ---")
    print(f"The time 't' to fall is given by the formula: t = sqrt(2 * h / g)")
    print(f"We use our calculated height h = {h} m and approximate gravity g = {g_approx} m/s^2.")
    
    t_squared_val = (2 * h) / g_approx
    print(f"First, we find t-squared: t^2 = (2 * {h}) / {g_approx} = {2*h} / {g_approx} = {int(t_squared_val)}")
    
    print(f"So, t = sqrt({int(t_squared_val)})")
    
    # --- Step 3: Simplify the square root ---
    # t_squared_val is 20. sqrt(20) = sqrt(4 * 5) = 2 * sqrt(5)
    print("To solve this, we simplify sqrt(20) to 2 * sqrt(5).")
    print(f"We approximate sqrt(5) with the fraction: {sqrt5_approx_num} / {sqrt5_approx_den}")
    
    final_time_num = 2 * sqrt5_approx_num
    final_time_den = sqrt5_approx_den
    
    # Simplify the final fraction 18/4 to 9/2
    common_divisor = math.gcd(final_time_num, final_time_den)
    simplified_num = final_time_num // common_divisor
    simplified_den = final_time_den // common_divisor

    print(f"Calculation: t ≈ 2 * ({sqrt5_approx_num} / {sqrt5_approx_den}) = {final_time_num}/{final_time_den} = {simplified_num}/{simplified_den} seconds.")
    print("\n")

    # --- Final Equation ---
    print("--- Final Equation Summary ---")
    print("The whole calculation in one line looks like this:")
    # Note: We print each number that goes into the final equation as requested.
    print(f"Time = sqrt( (2 * {distance_from_base} * {tan_alpha_num} / {tan_alpha_den}) / {g_approx} ) ≈ {simplified_num}/{simplified_den} seconds")

solve_egg_drop_problem()