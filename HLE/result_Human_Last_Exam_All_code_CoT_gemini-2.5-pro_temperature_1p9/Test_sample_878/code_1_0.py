import math

def solve_egg_drop_problem():
    """
    Calculates the time for an egg to fall from a skyscraper,
    presenting the steps using simple integer and fractional arithmetic.
    """
    print("Yes, it is possible to calculate the fall time. Here is a way to do it with simple numbers:")
    print("-" * 60)
    
    # --- Step 1: Calculate the height of the skyscraper ---
    print("Step 1: Find the height 'h' of the skyscraper.")
    print("The angle is a quarter of a right angle, so it is 90 / 4 = 22.5 degrees.")
    print("We can find the height using: h = distance * tan(angle)")
    print("For our angle, tan(22.5) is exactly sqrt(2) - 1.")
    print("To keep numbers simple, we can say sqrt(2) is very close to the fraction 17 / 12.")
    print("\nLet's calculate h:")
    print("h = 240 * (sqrt(2) - 1)")
    print("h ≈ 240 * (17 / 12 - 1)")
    print("  = 240 * (5 / 12)")
    print("  = (240 / 12) * 5")
    print("  = 20 * 5")
    h = 100
    print(f"  = {h} meters")
    
    print("-" * 60)
    
    # --- Step 2: Calculate the time squared (t^2) ---
    print("Step 2: Find the time 't' it takes for the egg to fall.")
    print("The formula is t^2 = (2 * h) / g, where 'g' is from gravity.")
    print("We will use h = 100 meters.")
    print("Gravity, g, is 9.8 m/s^2. As a fraction, this is 98 / 10, or 49 / 5.")
    print("\nLet's calculate t^2:")
    print("t^2 ≈ (2 * 100) / (49 / 5)")
    print("    = 200 / (49 / 5)")
    print("    = (200 * 5) / 49")
    numerator_t2 = 1000
    denominator_t2 = 49
    print(f"    = {numerator_t2} / {denominator_t2}")

    print("-" * 60)

    # --- Step 3: Find the final value for t ---
    print("Step 3: Calculate the final number for t by taking the square root.")
    print("t ≈ sqrt(1000 / 49)")
    print("  = sqrt(1000) / sqrt(49)")
    print("  = sqrt(100 * 10) / 7")
    print("  = (10 * sqrt(10)) / 7")
    print("\nNow we need a simple fraction for sqrt(10). A good one is 19 / 6.")
    print("t ≈ (10 / 7) * (19 / 6)")
    print("  = (10 * 19) / (7 * 6)")
    print("  = 190 / 42")
    print("We can simplify this fraction by dividing the top and bottom by 2:")
    final_numerator = 95
    final_denominator = 21
    print(f"t ≈ {final_numerator} / {final_denominator} seconds")

    # Verification of the error
    h_actual = 240 * (math.sqrt(2) - 1)
    g_actual = 9.8
    t_actual = math.sqrt(2 * h_actual / g_actual)
    t_approx = final_numerator / final_denominator
    error = abs(t_actual - t_approx)

    print("-" * 60)
    print(f"This answer is very close! The error is about {error:.2f} seconds, which is less than 0.1 seconds.")

solve_egg_drop_problem()