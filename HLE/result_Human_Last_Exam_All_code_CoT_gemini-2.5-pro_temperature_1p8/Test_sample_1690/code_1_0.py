import math

def solve_wuxing_problem():
    """
    Solves the Pioneer probe power problem based on the Wuxing architecture specifications.
    """
    # Given parameters
    P = 100.0  # Initial power (W)
    x = 4.0    # Control subsystem power (W)
    y = 5.0    # Camera subsystem power (W)
    z = 1.0    # Sensor subsystem power (W)
    T = 400.0  # Half-life in days

    # --- Part a: Calculate when to turn off the camera ---
    # This occurs when power drops below the requirement for all three subsystems.
    # Formula: t = T * (ln(P_initial) - ln(P_threshold)) / ln(2)
    power_threshold_a = x + y + z
    # t_a = T * math.log(P / power_threshold_a) / math.log(2)
    # Using the property log(a/b) = log(a) - log(b) to match the formula derivation
    t_a = T * (math.log(P) - math.log(power_threshold_a)) / math.log(2)
    
    # Round to the nearest total day
    answer_a = round(t_a)

    # --- Part b: Calculate when to turn off the sensor ---
    # This occurs when power drops below the requirement for control and sensor subsystems.
    power_threshold_b = x + z
    t_b = T * (math.log(P) - math.log(power_threshold_b)) / math.log(2)
    
    # Round to the nearest total day
    answer_b = round(t_b)
    
    # --- Part c: Calculate memory usage in D ---
    # In a memory-efficient program, variables would be:
    # 4 'int' variables for inputs (P, x, y, z): 4 * 5D = 20D
    # 1 'frac' variable for the constant T: 1 * 6D = 6D
    # 1 'frac' for pre-calculated ln(2): 1 * 6D = 6D
    # 1 'frac' for pre-calculated ln(P): 1 * 6D = 6D
    # 1 'frac' for the result, which can be reused: 1 * 6D = 6D
    # Total memory = 20D + 6D + 6D + 6D + 6D = 44D
    answer_c = 44

    # --- Part d: Calculate the number of calls to exp() ---
    # The program makes 4 ln() calls: ln(2), ln(P=100), ln(x+y+z=10), ln(x+z=5)
    # Each ln(val) call involves a binary search that calls exp() k times.
    # The number of iterations k is the smallest integer such that: val / 2^k < 1e-9
    # This can be rewritten as: k > log2(val * 1e9)
    # So, k = ceil(log2(val) + log2(1e9))
    
    log2_1e9 = math.log2(1e9)
    
    calls_for_ln_2 = math.ceil(math.log2(2) + log2_1e9)
    calls_for_ln_100 = math.ceil(math.log2(P) + log2_1e9)
    calls_for_ln_10 = math.ceil(math.log2(power_threshold_a) + log2_1e9)
    calls_for_ln_5 = math.ceil(math.log2(power_threshold_b) + log2_1e9)
    
    answer_d = calls_for_ln_2 + calls_for_ln_100 + calls_for_ln_10 + calls_for_ln_5

    # Print the final result in the format a:b:c:d
    print(f"{answer_a}:{answer_b}:{answer_c}:{answer_d}")

solve_wuxing_problem()
<<<1329:1729:44:135>>>