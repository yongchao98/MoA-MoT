import math

def solve_egg_drop_problem():
    """
    Calculates the egg drop time using simple integer and fraction arithmetic,
    and determines the largest integer used in the calculation.
    """
    # --- Setup and Approximations ---
    # Known distance from the skyscraper
    distance = 240

    # We are given that the angle alpha is 1/4 of a right angle, so tan(alpha) = tan(22.5) = sqrt(2) - 1.
    # We approximate tan(alpha) with a simple fraction to keep numbers small.
    # 5/12 is a good approximation for sqrt(2) - 1.
    tan_alpha_num = 5
    tan_alpha_den = 12

    # We use a simple integer approximation for gravity, g.
    g_approx = 10

    # The calculation will lead to sqrt(20). We approximate this with a fraction.
    # 9/2 = 4.5, and 4.5^2 = 20.25, which is very close to 20.
    sqrt_20_num = 9
    sqrt_20_den = 2

    # Keep track of the largest integer used in the calculation steps
    largest_integer = 0
    all_integers = set()

    # --- Step 1: Calculate the height 'h' ---
    print("Here is a simple way to calculate the time:")
    print("\nStep 1: Find the height of the skyscraper (h).")
    print(f"We approximate tan(alpha) with the fraction {tan_alpha_num}/{tan_alpha_den}.")
    
    # h = distance * tan(alpha)
    h_calc_result = distance * tan_alpha_num / tan_alpha_den
    h = int(h_calc_result)
    
    print(f"h = {distance} * {tan_alpha_num} / {tan_alpha_den} = {h}")
    
    # Update largest integer
    all_integers.update([distance, tan_alpha_num, tan_alpha_den, h])

    # --- Step 2: Calculate the square of the fall time (t^2) ---
    print("\nStep 2: Find the square of the time (t^2) using the formula t^2 = 2 * h / g.")
    print(f"We use g = {g_approx} for gravity.")
    
    # t^2 = 2 * h / g
    t_squared_calc_result = 2 * h / g_approx
    t_squared = int(t_squared_calc_result)
    
    print(f"t^2 = (2 * {h}) / {g_approx} = {t_squared}")
    
    # Update largest integer
    all_integers.update([2, h, g_approx, t_squared])

    # --- Step 3: Calculate the time 't' ---
    print(f"\nStep 3: Find the time (t) by taking the square root of {t_squared}.")
    print(f"We can approximate the square root of {t_squared} with the fraction {sqrt_20_num}/{sqrt_20_den}.")
    
    # t = sqrt(t_squared)
    final_time = sqrt_20_num / sqrt_20_den
    
    print(f"t = {sqrt_20_num} / {sqrt_20_den} = {final_time} seconds")
    
    # Update largest integer
    all_integers.update([sqrt_20_num, sqrt_20_den])
    
    largest_integer = max(all_integers)

    # --- Verification (for our reference) ---
    actual_h = 240 * (math.sqrt(2) - 1)
    actual_t = math.sqrt(2 * actual_h / 9.80665)
    error = abs(actual_t - final_time)
    
    print(f"\nThis calculation is possible for your son.")
    print(f"The final answer of {final_time}s has an error of about {error:.2f}s, which is less than the 0.1s requirement.")
    print(f"The largest integer that appears in this calculation is {largest_integer}.")

    # --- Final Answer ---
    print(f"\n<<<Y{largest_integer}>>>")

solve_egg_drop_problem()