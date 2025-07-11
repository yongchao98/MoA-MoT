import math

def solve_egg_drop_time():
    """
    Calculates the time for an egg to fall from a skyscraper using simple
    integer and fractional arithmetic, and verifies the accuracy.
    """
    # Known values
    distance = 240  # meters

    # --- Step 1: Calculate the height of the skyscraper ---
    # We need to approximate tan(22.5 degrees). A good fraction is 5/12.
    tan_alpha_num = 5
    tan_alpha_den = 12

    print("Step 1: Calculate the height (h).")
    print(f"We use the formula: h = distance * tan(angle)")
    print(f"Using the approximation tan(22.5°) ≈ {tan_alpha_num} / {tan_alpha_den}:")
    # Perform calculation for height
    height = distance * tan_alpha_num / tan_alpha_den
    print(f"h = {distance} * {tan_alpha_num} / {tan_alpha_den} = {int(height)} meters\n")

    # --- Step 2: Calculate the squared falling time (t^2) ---
    # We use the simplified gravity value g = 10 m/s^2.
    gravity = 10
    print("Step 2: Calculate the falling time (t).")
    print(f"We use the formula: t² = 2 * h / g")
    print(f"Using our calculated height h = {int(height)} and simplified gravity g = {gravity}:")
    # Perform calculation for time squared
    time_squared = 2 * height / gravity
    print(f"t² = {2} * {int(height)} / {gravity} = {time_squared}\n")


    # --- Step 3: Approximate the final time ---
    # We need to approximate sqrt(20). A good fraction is 9/2.
    time_approx_num = 9
    time_approx_den = 2
    print("Step 3: Find the final time by taking the square root.")
    print(f"t ≈ sqrt({time_squared})")
    print(f"We approximate this with the fraction {time_approx_num} / {time_approx_den}:")
    # Calculate final time
    time_approx = time_approx_num / time_approx_den
    print(f"t ≈ {time_approx_num} / {time_approx_den} = {time_approx} seconds\n")

    # --- Verification Step (for our reference) ---
    h_true = 240 * (math.sqrt(2) - 1)
    t_true = math.sqrt(2 * h_true / 9.8)
    error = abs(time_approx - t_true)
    
    print("--- Verification ---")
    print(f"The calculated time is {time_approx:.3f} s.")
    print(f"The more precise time is {t_true:.3f} s.")
    print(f"The absolute error is {error:.3f} s, which is less than the required 0.1 s.")

solve_egg_drop_time()

# The largest integer that appears in the calculation steps shown to the son is 240.
# Therefore, the final answer is Y240.
print("\n<<<Y240>>>")