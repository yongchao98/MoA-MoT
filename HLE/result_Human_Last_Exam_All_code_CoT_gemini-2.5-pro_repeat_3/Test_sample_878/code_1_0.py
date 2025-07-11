import math

def solve_egg_drop():
    """
    Calculates the egg drop time using simple integer and fraction arithmetic,
    explains the process, and determines the final answer code.
    """
    # Initial values from the problem
    d = 240 # distance from the skyscraper in meters

    # --- Approximations to keep numbers simple for a child ---
    # tan(22.5°) = sqrt(2) - 1 ≈ 0.414. A good fractional approximation is 5/12 ≈ 0.417.
    tan_alpha_num = 5
    tan_alpha_den = 12
    # Gravity g ≈ 9.8 m/s². A simple approximation is 10.
    g_approx = 10
    # For the final step, we will need to approximate sqrt(5). A good fraction is 9/4 = 2.25 (actual is ~2.236).
    sqrt5_approx_num = 9
    sqrt5_approx_den = 4

    print("Yes, your son can calculate the time! Here is a simple way to do it:")
    print("-" * 30)

    # --- Step 1: Calculate the height of the skyscraper ---
    h = d * tan_alpha_num / tan_alpha_den
    print("Step 1: Find the height (h) of the skyscraper.")
    print(f"The distance is {d}m. We approximate the tangent of the angle as {tan_alpha_num}/{tan_alpha_den}.")
    print(f"The calculation for the height is:")
    print(f"h = {d} * {tan_alpha_num} / {tan_alpha_den} = {int(h)} meters")
    print("-" * 30)

    # --- Step 2: Use the height to find the time ---
    t_squared = (2 * h) / g_approx
    print("Step 2: Use the height to find the fall time (t).")
    print(f"The formula is h = 1/2 * g * t^2. We can use g = {g_approx} for simplicity.")
    print(f"Rearranging for t^2 gives: t^2 = 2 * h / g")
    print(f"The calculation for t^2 is:")
    print(f"t^2 = 2 * {int(h)} / {g_approx} = {int(t_squared)}")
    print("-" * 30)

    # --- Step 3: Simplify the square root ---
    # t = sqrt(20) = sqrt(4*5) = 2 * sqrt(5)
    t_final_num = 2 * sqrt5_approx_num
    t_final_den = sqrt5_approx_den
    t_final_val = t_final_num / t_final_den

    print("Step 3: Calculate the final time by taking the square root.")
    print(f"So, t = sqrt({int(t_squared)}).")
    print(f"We can simplify sqrt({int(t_squared)}) to 2 * sqrt(5).")
    print(f"A good fraction for sqrt(5) is {sqrt5_approx_num}/{sqrt5_approx_den}.")
    print(f"The final calculation for time is:")
    print(f"t = 2 * {sqrt5_approx_num} / {sqrt5_approx_den} = {t_final_num / t_final_den} seconds")
    print("-" * 30)

    # --- Final Answer Generation ---
    # The integers used in the printed calculation are:
    # 240, 5, 12, 100, 2, 10, 20, 9, 4
    # The largest integer is 240.
    # The problem asks for Yz where z is the largest integer.
    largest_integer = 240
    answer_code = f"Y{largest_integer}"

    # This part is just to verify the error is less than 0.1s
    # t_exact = math.sqrt(2 * 240 * (math.sqrt(2) - 1) / 9.8)
    # error = abs(t_final_val - t_exact)
    # print(f"(Self-check: Exact time is ~{t_exact:.3f}s, error is ~{error:.3f}s, which is < 0.1s)")
    
    print(f"<<<{answer_code}>>>")

solve_egg_drop()