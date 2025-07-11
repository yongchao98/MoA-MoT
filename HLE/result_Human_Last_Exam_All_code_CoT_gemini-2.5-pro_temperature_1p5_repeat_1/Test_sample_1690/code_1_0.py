import math

def solve_pioneer_problem():
    """
    Calculates the answers to the Pioneer probe problem.
    """
    # Given parameters
    P = 100  # Initial power in W
    x = 4    # Control subsystem power in W
    y = 5    # Camera subsystem power in W
    z = 1    # Sensor subsystem power in W
    half_life = 400  # Half-life of the battery in days

    # --- Part a: Camera subsystem shutdown time ---
    # The power threshold is when power is needed for all systems.
    threshold_a = x + y + z
    # The time is calculated using the formula: t = T * ln(P / P_threshold) / ln(2)
    time_a = half_life * math.log(P / threshold_a) / math.log(2)
    # The result is rounded to the nearest whole day.
    a = round(time_a)

    print(f"a. The camera is turned off when power drops below {threshold_a} W.")
    print(f"   Equation: t = {half_life} * (ln({P}) - ln({threshold_a})) / ln(2)")
    print(f"   The numbers in the equation are: {half_life}, {P}, {threshold_a}, and 2.")
    print(f"   Result: {a} days\n")

    # --- Part b: Sensor subsystem shutdown time ---
    # The power threshold is when power is needed for control and sensor systems.
    threshold_b = x + z
    # The time is calculated using the same formula with the new threshold.
    time_b = half_life * math.log(P / threshold_b) / math.log(2)
    # The result is rounded to the nearest whole day.
    b = round(time_b)
    
    print(f"b. The sensor is turned off when power drops below {threshold_b} W.")
    print(f"   Equation: t = {half_life} * (ln({P}) - ln({threshold_b})) / ln(2)")
    print(f"   The numbers in the equation are: {half_life}, {P}, {threshold_b}, and 2.")
    print(f"   Result: {b} days\n")

    # --- Part c: Memory usage for the most efficient Wuxing C program ---
    # The most time-efficient program minimizes calls to ln() by computing ln(2) and ln(5)
    # and deriving other values.
    # Variables required:
    # - frac ln2 (6D)
    # - frac ln5 (6D)
    # - frac t (reusable for t_a and t_b) (6D)
    # - int days_a (5D)
    # - int days_b (5D)
    # Total memory = 3 * 6D (fracs) + 2 * 5D (ints) = 18D + 10D = 28D.
    c = 28
    
    print(f"c. The memory usage of the most efficient program is {c} D.\n")

    # --- Part d: Number of calls to the exp() function ---
    # The efficient program only calls ln(2) and ln(5).
    # The number of exp() calls within ln(x) is ceil(log2(x * 10^9)).
    exp_calls_for_ln2 = math.ceil(math.log2(2 * 10**9))
    exp_calls_for_ln5 = math.ceil(math.log2(5 * 10**9))
    d = exp_calls_for_ln2 + exp_calls_for_ln5
    
    print(f"d. The number of calls to the expensive exp() function is {d}.\n")
    
    # --- Final formatted answer ---
    final_answer = f"{a}:{b}:{c}:{d}"
    print(f"Final Answer: {final_answer}")
    return final_answer

# Execute the function
solve_pioneer_problem()