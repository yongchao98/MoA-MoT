import math

def solve_pioneer_problem():
    """
    Calculates the answers for the Pioneer probe problem based on Wuxing architecture constraints.
    """
    # Part 1 & 2: Define parameters and formula
    P = 100.0  # Initial power (W)
    x = 4.0    # Control subsystem power (W)
    y = 5.0    # Camera subsystem power (W)
    z = 1.0    # Sensor subsystem power (W)
    T_half = 400.0  # Half-life (days)

    # The core equation to find time 't' is:
    # t = T_half * ln(P_target / P_initial) / ln(0.5)

    # --- Part a: Camera subsystem shutdown ---
    print("--- Calculating 'a': Camera shutdown time ---")
    power_req_a = x + y + z
    ratio_a = power_req_a / P
    t_a_float = T_half * math.log(ratio_a) / math.log(0.5)
    t_a_rounded = round(t_a_float)
    
    print(f"Power threshold (x+y+z): {x} + {y} + {z} = {power_req_a} W")
    print(f"Equation: t = {T_half} * ln({power_req_a} / {P}) / ln(0.5)")
    print(f"Result: t = {t_a_float:.2f} days, rounded to {t_a_rounded} days.")
    print(f"a = {t_a_rounded}\n")


    # --- Part b: Sensor subsystem shutdown ---
    print("--- Calculating 'b': Sensor shutdown time ---")
    power_req_b = x + z
    ratio_b = power_req_b / P
    t_b_float = T_half * math.log(ratio_b) / math.log(0.5)
    t_b_rounded = round(t_b_float)
    
    print(f"Power threshold (x+z): {x} + {z} = {power_req_b} W")
    print(f"Equation: t = {T_half} * ln({power_req_b} / {P}) / ln(0.5)")
    print(f"Result: t = {t_b_float:.2f} days, rounded to {t_b_rounded} days.")
    print(f"b = {t_b_rounded}\n")


    # --- Part c: Memory usage analysis ---
    print("--- Calculating 'c': Memory usage in D ---")
    # An optimized program would declare variables for:
    # Inputs (P,x,y,z), Constants (T_half, half), Intermediates (ln_half, temp_frac) -> 8 frac types
    # Results (t_a, t_b) -> 2 int types
    frac_vars = 8
    int_vars = 2
    frac_size_d = 6 # 2D(n) + 2D(d) + 2D(e)
    int_size_d = 5
    memory_usage_c = (frac_vars * frac_size_d) + (int_vars * int_size_d)

    print(f"Optimized variables: {frac_vars} of type 'frac', {int_vars} of type 'int'")
    print(f"Memory calculation: ({frac_vars} * {frac_size_d}D) + ({int_vars} * {int_size_d}D)")
    print(f"Total memory usage: {memory_usage_c}D")
    print(f"c = {memory_usage_c}\n")


    # --- Part d: exp function call analysis ---
    print("--- Calculating 'd': Number of 'exp' calls ---")
    # Total calls is the sum of calls from ln(0.5), ln(ratio_a), and ln(ratio_b)
    # Number of calls for ln(arg) is k = ceil(log2(arg * 10^9))
    calls_ln_half = math.ceil(math.log2(0.5 * 10**9))
    calls_ln_a = math.ceil(math.log2(ratio_a * 10**9))
    calls_ln_b = math.ceil(math.log2(ratio_b * 10**9))
    total_exp_calls = calls_ln_half + calls_ln_a + calls_ln_b
    
    print(f"Calls for ln(0.5): ceil(log2(0.5e9)) = {calls_ln_half}")
    print(f"Calls for ln({ratio_a}): ceil(log2({ratio_a}e9)) = {calls_ln_a}")
    print(f"Calls for ln({ratio_b}): ceil(log2({ratio_b}e9)) = {calls_ln_b}")
    print(f"Total 'exp' calls: {calls_ln_half} + {calls_ln_a} + {calls_ln_b} = {total_exp_calls}")
    print(f"d = {total_exp_calls}\n")
    
    # --- Final result in specified format ---
    print(f"Final Answer (a:b:c:d): {t_a_rounded}:{t_b_rounded}:{memory_usage_c}:{total_exp_calls}")

solve_pioneer_problem()