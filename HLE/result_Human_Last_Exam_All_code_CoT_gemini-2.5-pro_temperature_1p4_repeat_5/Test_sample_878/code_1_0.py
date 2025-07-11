import math

def solve_for_son():
    """
    Calculates the fall time of an egg from a skyscraper using simple integer
    and fraction arithmetic suitable for a child.
    """

    print("Yes, your son can solve this! Here is a way to do the calculation with simple numbers.")
    print("-" * 50)
    
    # --- STEP 1: Calculate the height of the skyscraper ---
    print("Step 1: Let's find the height of the skyscraper (h).")
    
    d = 240  # distance from the skyscraper in meters
    
    # The angle alpha is 1/4 of a right angle (90 degrees), so alpha = 22.5 degrees.
    # The exact value of tan(22.5) is sqrt(2) - 1.
    # To keep numbers simple, we'll use a fraction for tan(alpha).
    # A good fraction for tan(22.5) is 5/12.
    tan_alpha_num = 5
    tan_alpha_den = 12
    
    # Calculate height h
    h = d * tan_alpha_num / tan_alpha_den
    
    print(f"The height 'h' is the distance '{d} m' multiplied by a special number for the angle.")
    print(f"For our angle, we can use the simple fraction: {tan_alpha_num}/{tan_alpha_den}")
    print("\nCalculation for height:")
    print(f"h = {d} * {tan_alpha_num} / {tan_alpha_den}")
    print(f"h = {d * tan_alpha_num} / {tan_alpha_den}")
    print(f"h = {int(h)} meters")
    
    print("-" * 50)
    
    # --- STEP 2: Calculate the time to fall ---
    print("Step 2: Now let's find the time (t) it takes for the egg to fall.")
    
    # Use a simple approximation for gravity, g = 10 m/s^2
    g = 10
    
    # From h = (1/2)*g*t^2, we get t^2 = 2*h/g
    t_squared = 2 * h / g
    
    print(f"We use the rule: t * t = 2 * h / g")
    print(f"For gravity 'g', we can use the simple number {g}.")
    print("\nCalculation for time squared (t*t):")
    print(f"t * t = (2 * {int(h)}) / {g}")
    print(f"t * t = {2 * int(h)} / {g}")
    print(f"t * t = {int(t_squared)}")
    
    # Now we need to find the square root of t_squared (20)
    # A good fractional approximation for sqrt(20) is 9/2, since (9/2)*(9/2) = 81/4 = 20.25
    t_approx_num = 9
    t_approx_den = 2
    
    print(f"\nTo find 't', we need a number that, when multiplied by itself, is close to {int(t_squared)}.")
    print(f"Let's try the fraction {t_approx_num}/{t_approx_den}. See: {t_approx_num}/{t_approx_den} * {t_approx_num}/{t_approx_den} = {t_approx_num**2}/{t_approx_den**2}, which is very close to {int(t_squared)}!")
    print(f"So, the time 't' is about {t_approx_num / t_approx_den} seconds.")
    
    print("-" * 50)

    # --- STEP 3: Final Combined Equation ---
    print("Here is the whole calculation put together in one line:")
    print(f"Time ≈ sqrt( (2 * ({d} * {tan_alpha_num} / {tan_alpha_den})) / {g} ) ≈ {t_approx_num}/{t_approx_den}")
    
    # --- Error Check (for verification, not shown to son) ---
    h_actual = 240 * (math.sqrt(2) - 1)
    t_actual = math.sqrt(2 * h_actual / 9.8)
    t_approx = t_approx_num / t_approx_den
    error = abs(t_actual - t_approx)
    # print(f"\n[Verification: The error is {error:.4f}s, which is less than 0.1s]")


solve_for_son()
<<<Y240>>>