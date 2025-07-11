import math

def solve_pioneer_problem():
    """
    Solves the Pioneer probe power problem based on the Wuxing architecture.
    """

    # --- Given values ---
    P_initial = 100.0  # Initial power in Watts
    x = 4.0            # Control subsystem power
    y = 5.0            # Camera subsystem power
    z = 1.0            # Sensor subsystem power
    T_half = 400.0     # Power half-life in days

    # --- Part a: Camera subsystem shutdown time ---
    # Shutdown occurs when power P(t) drops below the total required power.
    p_threshold_a = x + y + z
    # The formula to find time t is: t = T_half * ln(P_initial / P_threshold) / ln(2)
    t_a_raw = T_half * math.log(P_initial / p_threshold_a) / math.log(2)
    a = round(t_a_raw)

    # --- Part b: Sensor subsystem shutdown time ---
    # Shutdown occurs when power P(t) drops below the power for control and sensor.
    p_threshold_b = x + z
    t_b_raw = T_half * math.log(P_initial / p_threshold_b) / math.log(2)
    b = round(t_b_raw)
    
    # --- Part c: Memory usage for the most efficient program ---
    # The most time-efficient program uses the relation t_b = T_half + t_a.
    # It requires ln(2) and ln(P/p_threshold_a) = ln(10).
    # Variables needed: P, x, y, z (inputs, 4 frac), T_half (constant, 1 frac),
    # P_div_Pa (intermediate, 1 frac), ln_2 (1 frac), ln_P_div_Pa (1 frac),
    # t_a (result, 1 frac), t_b (result, 1 frac).
    # Total frac variables = 4 + 1 + 1 + 1 + 1 + 1 + 1 = 10
    num_frac_variables = 10
    # Each char is 2D, frac has 3 chars. So, size of frac = 3 * 2D = 6D.
    size_of_frac_in_D = 6
    c = num_frac_variables * size_of_frac_in_D

    # --- Part d: Number of 'exp' calls in the most efficient program ---
    # The efficient program only makes two ln calls: ln(2) and ln(10).
    # The number of calls to exp() inside ln(val) is ceil(log2(val / 1e-9)).
    
    # Number of calls for ln(2)
    ln_arg_2 = 2.0
    exp_calls_for_ln2 = math.ceil(math.log2(ln_arg_2 / 1e-9))
    
    # Number of calls for ln(10)
    ln_arg_10 = 10.0
    exp_calls_for_ln10 = math.ceil(math.log2(ln_arg_10 / 1e-9))
    
    d = exp_calls_for_ln2 + exp_calls_for_ln10

    # Print the final result in the format a:b:c:d
    # The prompt also asks to output each number in the final equation.
    print(f"a (camera shutdown day): {a}")
    print(f"b (sensor shutdown day): {b}")
    print(f"c (memory usage in D): {c}")
    print(f"d (number of exp calls): {d}")
    print(f"\nFinal Answer String:")
    print(f"{a}:{b}:{c}:{d}")


solve_pioneer_problem()