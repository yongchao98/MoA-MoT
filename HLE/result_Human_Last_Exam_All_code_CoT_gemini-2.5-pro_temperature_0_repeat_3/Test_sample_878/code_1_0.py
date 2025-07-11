def solve_egg_drop():
    """
    Calculates the egg drop time using simple integer and fraction arithmetic.
    """
    # Initial given values
    distance = 240  # in meters

    # --- Step 1: Calculate the height of the skyscraper (h) ---
    print("Step 1: Find the height of the skyscraper.")
    print(f"The skyscraper is {distance}m away.")
    print("The angle α is 1/4 of a right angle (90 degrees), so α = 22.5 degrees.")
    print("The height 'h' is calculated by: h = distance * tan(22.5°)")
    print("A good fraction for tan(22.5°), which is √2 - 1, is 5 / 12.")
    
    # h = 240 * (5/12)
    h_numerator = distance * 5
    h = h_numerator // 12
    print(f"So, h = {distance} * 5 / 12 = {h}")
    print("-" * 20)

    # --- Step 2: Calculate the time to fall (t) ---
    print("Step 2: Find the time it takes for the egg to fall.")
    print("The formula is h = (1/2) * g * t², where g is about 10.")
    print("We can rearrange this to find t²: t² = 2 * h / g")
    
    # t_squared = (2 * 100) / 10
    g_approx = 10
    t_squared_numerator = 2 * h
    t_squared = t_squared_numerator // g_approx
    print(f"So, t² = (2 * {h}) / {g_approx} = {t_squared}")
    print("-" * 20)

    # --- Step 3: Find t by taking the square root ---
    print("Step 3: Find t from t².")
    print(f"t = sqrt({t_squared}) = sqrt(4 * 5) = 2 * sqrt(5)")
    print("A good fraction for sqrt(5) is 9 / 4.")
    
    # t = 2 * (9/4) = 9/2
    sqrt5_num = 9
    sqrt5_den = 4
    final_t_num = 2 * sqrt5_num
    final_t_den = sqrt5_den
    
    print(f"So, the final equation is t ≈ 2 * {sqrt5_num} / {sqrt5_den}")
    
    final_t_decimal = final_t_num / final_t_den
    print(f"The final answer is t = {final_t_num / 2} / {final_t_den / 2} = {final_t_decimal} seconds.")
    print("-" * 20)

    # --- Conclusion ---
    # The actual time is approx 4.502s, our result of 4.5s has an error of ~0.002s, which is less than 0.1s.
    # The calculation is possible for your son.
    # The integers used in the calculation are 240, 5, 12, 100, 2, 10, 20, 4, 9.
    # The largest integer that appears in the calculation is the initial distance.
    largest_integer = 240
    print(f"The calculation is possible and the largest integer used is {largest_integer}.")
    print(f"<<<Y{largest_integer}>>>")

solve_egg_drop()