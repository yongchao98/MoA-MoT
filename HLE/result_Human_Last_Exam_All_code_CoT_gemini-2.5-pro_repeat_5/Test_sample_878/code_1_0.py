def solve_egg_drop_problem():
    """
    This function calculates the egg drop time using simple integer and
    fractional arithmetic, as requested.
    """

    # --- Introduction for the young mathematician ---
    print("Let's figure out how long it takes for the egg to fall!")
    print("We will solve this in three simple steps.\n")

    # --- Step 1: Calculate the height of the skyscraper (h) ---
    print("--- Step 1: Find the height of the skyscraper ---")
    d = 240  # distance from the skyscraper in meters
    # We approximate tan(22.5 degrees) as 5/12
    tan_alpha_num = 5
    tan_alpha_den = 12
    h = d * tan_alpha_num / tan_alpha_den
    print(f"The skyscraper is {d}m away.")
    print(f"We can find the height 'h' with the equation: h = {d} * ({tan_alpha_num}/{tan_alpha_den})")
    print(f"So, the height calculation is:")
    print(f"h = {d} * {tan_alpha_num} / {tan_alpha_den} = {int(h)} meters")
    print("-" * 20 + "\n")

    # --- Step 2: Set up the equation for the fall time (t) ---
    print("--- Step 2: Find the squared fall time (t^2) ---")
    # We approximate gravity g as 10 m/s^2
    g = 10
    # The formula is t^2 = 2 * h / g
    t_squared = 2 * h / g
    print(f"The formula for the fall time 't' is t^2 = (2 * h) / g.")
    print(f"We use our height h = {int(h)} and approximate gravity g = {g}.")
    print(f"So, the calculation for t-squared is:")
    print(f"t^2 = (2 * {int(h)}) / {g} = {int(t_squared)}")
    print("-" * 20 + "\n")

    # --- Step 3: Approximate the final time (t) ---
    print("--- Step 3: Find the final time 't' ---")
    # We need to find the square root of t_squared, which is 20.
    # We approximate sqrt(20) as 9/2, because (9/2)^2 = 81/4 = 20.25, which is very close to 20.
    final_t_num = 9
    final_t_den = 2
    final_t = final_t_num / final_t_den
    print(f"Now we just need to find the square root of {int(t_squared)}.")
    print(f"A good and simple fraction for this is {final_t_num}/{final_t_den}.")
    print(f"So, the final calculation for the time 't' is:")
    print(f"t ≈ {final_t_num} / {final_t_den} = {final_t} seconds")

solve_egg_drop_problem()
<<<Y240>>>