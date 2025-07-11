import math

def solve_egg_drop():
    """
    Prints a step-by-step guide to calculating the egg drop time
    using approximations with small integers and fractions.
    """
    # --- Introduction for the user's son ---
    print("Yes, you can definitely calculate this! Let's solve it together step by step, using small numbers.")
    print("----------------------------------------------------------------------------------------------------")
    print("Our goal is to find the time 't' an egg takes to fall from a skyscraper.\n")

    # --- Step 1: Calculate the height (h) ---
    print("Step 1: Find the height of the skyscraper.")
    distance = 240
    angle_divisor = 4
    print(f"We are standing {distance}m away from the skyscraper.")
    print(f"The angle to the top is one-fourth of a right angle (90 degrees), so the angle is 90 / {angle_divisor} = 22.5 degrees.")
    print("The height 'h' is given by the formula: h = distance * tan(angle).")
    print(f"So, h = {distance} * tan(22.5°).")
    print("\nTo make it easy, we can replace tan(22.5°) with a simple fraction. Its value is about 0.414, which is very close to 2/5.")
    h_num, h_den = 2, 5
    approx_h = distance * h_num / h_den
    print(f"Let's use this simple fraction: h ≈ {distance} * {h_num} / {h_den}")
    print(f"So, the height is approximately {int(approx_h)} meters.\n")

    # --- Step 2: Calculate the fall time (t) ---
    print("Step 2: Find the time 't' it takes to fall.")
    h = int(approx_h)
    g_num, g_den = 49, 5
    print("The physics formula for the time is t = square_root(2 * h / g).")
    print(f"For gravity 'g', we will use the fraction {g_num}/{g_den}, which is equal to 9.8.")
    print("\nLet's calculate the part inside the square root first: 2 * h / g")
    print(f"= (2 * {h}) / ({g_num}/{g_den})")
    intermediate_val = 2 * h
    print(f"= {intermediate_val} / ({g_num}/{g_den})")
    print(f"= ({intermediate_val} * {g_den}) / {g_num}")
    t_squared_num = intermediate_val * g_den
    t_squared_den = g_num
    print(f"= {t_squared_num} / {t_squared_den}\n")

    # --- Step 3: Simplify and get the final answer ---
    print("Step 3: Find the square root to get the final answer.")
    print(f"Now we need to calculate the square root of {t_squared_num} / {t_squared_den}.")
    print(f"Taking the square root of {t_squared_num} is tricky because it's not a perfect square.")
    print("However, look! 31 * 31 = 961, which is very close to 960.")
    approx_t_squared_num = 961
    print(f"Let's make our calculation easier by approximating {t_squared_num} with {approx_t_squared_num}.")
    
    final_t_num = int(math.sqrt(approx_t_squared_num))
    final_t_den = int(math.sqrt(t_squared_den))
    
    print(f"\nSo, our final equation is:")
    print(f"t ≈ square_root({approx_t_squared_num} / {t_squared_den}) = square_root({approx_t_squared_num}) / square_root({t_squared_den}) = {final_t_num} / {final_t_den}")
    print(f"\nThe time it takes for the egg to reach the ground is about {final_t_num}/{final_t_den} seconds!")
    
    # --- Error Check (for the parent) ---
    t_exact = math.sqrt(2 * 240 * (math.sqrt(2)-1) / 9.8)
    t_approx = final_t_num / final_t_den
    error = abs(t_exact - t_approx)
    print(f"\n(Note for dad: The calculated time is {t_approx:.2f}s, and the true time is about {t_exact:.2f}s. The error of {error:.2f}s is less than 0.1s, so this simple calculation works perfectly!)")

solve_egg_drop()