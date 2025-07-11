import math

def solve_and_explain():
    """
    This function explains and calculates the egg drop time using simple fractions.
    """
    print("Yes, your son can calculate the time with an error less than 0.1s.")
    print("Here is a step-by-step calculation designed for him, using small integers and fractions.")
    
    print("\n--- Step 1: Find the height (h) of the skyscraper ---")
    d = 240
    print(f"The distance to the skyscraper is d = {d} m.")
    print("The angle is a fourth of a right angle, so it's 90 / 4 = 22.5 degrees.")
    print("To find the height, we use the formula: h = d * tan(22.5°).")
    print("A clever way to write tan(22.5°) is sqrt(2) - 1.")
    print("We can use a simple fraction to approximate sqrt(2), which is almost 7 / 5.")
    sqrt2_num, sqrt2_den = 7, 5
    print(f"So, tan(22.5°) ≈ {sqrt2_num}/{sqrt2_den} - 1 = {sqrt2_num - sqrt2_den}/{sqrt2_den}.")
    tan_a_num, tan_a_den = 2, 5
    
    print("\nNow, let's calculate the height h:")
    h_val_num = d * tan_a_num
    h_val_den = tan_a_den
    h_val = h_val_num // h_val_den
    print(f"h ≈ {d} * ({tan_a_num} / {tan_a_den})")
    print(f"h ≈ {h_val_num} / {h_val_den} = {h_val} meters")

    print("\n--- Step 2: Find the time (t) for the egg to fall ---")
    print("The formula for a falling object is t^2 = (2 * h) / g.")
    print("We will use g ≈ 9.8 m/s², which is the fraction 98/10 or 49/5.")
    g_num, g_den = 49, 5
    
    print("\nLet's calculate t-squared (t^2):")
    print(f"t^2 ≈ (2 * {h_val}) / ({g_num} / {g_den})")
    # This is equivalent to (2 * h) * (g_den / g_num)
    t_sq_num = 2 * h_val * g_den
    t_sq_den = g_num
    print(f"t^2 ≈ (2 * {h_val} * {g_den}) / {g_num}")
    print(f"t^2 ≈ {t_sq_num} / {t_sq_den}")

    print("\n--- Step 3: Find the final answer for t ---")
    print(f"The number {t_sq_num} is very close to 961, which is a special number because 31 * 31 = 961.")
    t_sq_num_approx = 961
    print(f"So, we can say that t^2 is approximately {t_sq_num_approx} / {t_sq_den}.")
    
    t_num = 31
    t_den = 7
    print("Now we just need to take the square root of the numerator and the denominator.")
    print(f"t ≈ sqrt({t_sq_num_approx}) / sqrt({t_sq_den})")
    print("\nThe final calculation is:")
    print(f"t ≈ {t_num} / {t_den} seconds")

    # This part is to verify the error is less than 0.1s
    # t_approx = t_num / t_den
    # h_true = 240 * (math.sqrt(2) - 1)
    # g_true = 9.8
    # t_true = math.sqrt(2 * h_true / g_true)
    # error = abs(t_true - t_approx)
    # print(f"\nFor verification: Approximate time is {t_approx:.4f}s, Real time is {t_true:.4f}s. Error is {error:.4f}s.")

solve_and_explain()