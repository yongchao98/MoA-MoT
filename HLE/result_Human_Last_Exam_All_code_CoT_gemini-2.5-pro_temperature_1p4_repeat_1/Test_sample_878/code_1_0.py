import math

def solve_egg_drop_problem():
    """
    Calculates the egg drop time using simple integer and fractional approximations,
    and presents the calculation in a way a child can understand.
    """
    # --- Part 1: Define initial values and approximations ---
    # Given distance from the skyscraper
    distance = 240  # meters

    # Approximations for the calculation
    # tan(22.5°) = sqrt(2) - 1 ≈ 0.414. The fraction 5/12 = 0.4166... is a good approximation.
    tan_alpha_num = 5
    tan_alpha_den = 12

    # A common and simple approximation for gravity (g)
    g_approx = 10  # m/s^2

    # We will need to approximate sqrt(5) later in the calculation.
    # sqrt(5) ≈ 2.236. The fraction 9/4 = 2.25 is a good approximation.
    sqrt5_num = 9
    sqrt5_den = 4

    # --- Part 2: Perform the calculation with approximations ---
    # Step 1: Calculate the approximate height (h)
    h_approx = distance * tan_alpha_num / tan_alpha_den
    
    # Step 2: Calculate the approximate time squared (t^2)
    # The formula is t^2 = 2 * h / g
    t_squared_approx = (2 * h_approx) / g_approx

    # Step 3: Calculate the approximate time (t)
    # t = sqrt(t_squared_approx) = sqrt(20) = 2 * sqrt(5)
    # We use our approximation for sqrt(5)
    final_t_num = 2 * sqrt5_num
    final_t_den = sqrt5_den
    t_approx_val = final_t_num / final_t_den

    # --- Part 3: Verify the accuracy ---
    # Calculate the "true" time for comparison
    g_true = 9.8
    h_true = distance * (math.sqrt(2) - 1)
    t_true = math.sqrt(2 * h_true / g_true)
    error = abs(t_true - t_approx_val)

    # --- Part 4: Print the results ---
    if error < 0.1:
        print("Yes, the calculation is possible with an absolute error less than 0.1s.")
        print("\nHere is a calculation your son can follow:")
        
        print("\n1. First, we find the height of the skyscraper (h).")
        print(f"We can approximate tan(22.5 degrees) with the fraction {tan_alpha_num}/{tan_alpha_den}.")
        print(f"h = {distance} * {tan_alpha_num} / {tan_alpha_den} = {int(h_approx)} meters")
        
        print("\n2. Next, we find the time it takes to fall.")
        print(f"The formula for time squared is t^2 = 2 * h / g.")
        print(f"We can use a simple value for gravity, g = {g_approx}.")
        print(f"t^2 = 2 * {int(h_approx)} / {g_approx} = {int(t_squared_approx)}")

        print("\n3. Now, we find the time (t) by taking the square root.")
        print(f"t = square root of {int(t_squared_approx)}, which simplifies to 2 * square root of 5.")
        print(f"We can approximate the square root of 5 with the fraction {sqrt5_num}/{sqrt5_den}.")
        
        print("\nSo, the final calculation for time is:")
        print(f"t = 2 * {sqrt5_num} / {sqrt5_den} = {final_t_num} / {final_t_den}")

        # Determine the largest integer used in the calculation steps shown above
        integers_used = [distance, tan_alpha_num, tan_alpha_den, int(h_approx), 2, g_approx, int(t_squared_approx), 5, sqrt5_num, sqrt5_den, final_t_num, final_t_den]
        largest_integer = max(integers_used)
        print(f"\n<<<Y{largest_integer}>>>")
    else:
        print("No, this calculation cannot be done with the required accuracy using these simple numbers.")
        print("\n<<<N0>>>")

solve_egg_drop_problem()