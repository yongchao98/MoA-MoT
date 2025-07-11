import math

def solve_pioneer_problem():
    """
    Calculates and explains the solution to the Pioneer probe problem.
    """
    # Given values
    P = 100  # Initial power in Watts
    x = 4    # Control subsystem power in Watts
    y = 5    # Camera subsystem power in Watts
    z = 1    # Sensor subsystem power in Watts
    HALF_LIFE = 400 # days
    PRECISION = 1e-9

    # --- Part a: Camera subsystem shutdown time ---
    power_req_a = x + y + z
    time_a_unrounded = HALF_LIFE * math.log(P / power_req_a) / math.log(2)
    result_a = round(time_a_unrounded)

    print("a. Time to turn off the camera subsystem (in days)")
    print(f"   Equation: t = {HALF_LIFE} * ln(P / (x + y + z)) / ln(2)")
    print(f"   Calculation: t = {HALF_LIFE} * ln({P} / ({x} + {y} + {z})) / ln(2)")
    print(f"   Calculation: t = {HALF_LIFE} * ln({P / power_req_a}) / ln(2)")
    print(f"   Result: t = {time_a_unrounded:.2f} days")
    print(f"   Rounded to the nearest day (a): {result_a}")
    print("-" * 20)

    # --- Part b: Sensor subsystem shutdown time ---
    power_req_b = x + z
    time_b_unrounded = HALF_LIFE * math.log(P / power_req_b) / math.log(2)
    result_b = round(time_b_unrounded)

    print("b. Time to turn off the sensor subsystem (in days)")
    print(f"   Equation: t = {HALF_LIFE} * ln(P / (x + z)) / ln(2)")
    print(f"   Calculation: t = {HALF_LIFE} * ln({P} / ({x} + {z})) / ln(2)")
    print(f"   Calculation: t = {HALF_LIFE} * ln({P / power_req_b}) / ln(2)")
    print(f"   Result: t = {time_b_unrounded:.2f} days")
    print(f"   Rounded to the nearest day (b): {result_b}")
    print("-" * 20)

    # --- Part c: Memory Usage ---
    # Memory usage of the most efficient C program variables:
    # 4 int inputs (P,x,y,z): 4 * 5D = 20D
    # 2 int results (days_a, days_b): 2 * 5D = 10D
    # 4 frac variables (t, ratio, ln_2, ln_ratio): 4 * 6D = 24D
    D_INT = 5
    D_FRAC = 6
    num_input_ints = 4
    num_result_ints = 2
    num_frac_vars = 4
    result_c = (num_input_ints * D_INT) + (num_result_ints * D_INT) + (num_frac_vars * D_FRAC)

    print("c. Memory usage in D (decimal digits) for variables")
    print(f"   int P, x, y, z;      // {num_input_ints} * {D_INT}D = {num_input_ints * D_INT}D")
    print(f"   int days_a, days_b;  // {num_result_ints} * {D_INT}D = {num_result_ints * D_INT}D")
    print(f"   frac t, r, l2, lr; // {num_frac_vars} * {D_FRAC}D = {num_frac_vars * D_FRAC}D")
    print(f"   Total memory usage (c): {num_input_ints*D_INT} + {num_result_ints*D_INT} + {num_frac_vars*D_FRAC} = {result_c}D")
    print("-" * 20)

    # --- Part d: Number of 'exp' calls ---
    # The most efficient program makes 3 ln() calls: ln(2), ln(10), ln(20)
    # Number of exp calls for ln(val) is ceil(log2(val / precision))
    val_ln_1 = 2
    val_ln_2 = P / power_req_a # = 10
    val_ln_3 = P / power_req_b # = 20

    exp_calls_1 = math.ceil(math.log2(val_ln_1 / PRECISION))
    exp_calls_2 = math.ceil(math.log2(val_ln_2 / PRECISION))
    exp_calls_3 = math.ceil(math.log2(val_ln_3 / PRECISION))
    result_d = exp_calls_1 + exp_calls_2 + exp_calls_3
    
    print("d. Number of times the program calls function 'exp'")
    print(f"   The efficient program makes 3 ln() calls: ln({val_ln_1}), ln({val_ln_2}), ln({val_ln_3})")
    print(f"   exp_calls for ln({val_ln_1}) = ceil(log2({val_ln_1} / {PRECISION})) = {exp_calls_1}")
    print(f"   exp_calls for ln({val_ln_2}) = ceil(log2({val_ln_2} / {PRECISION})) = {exp_calls_2}")
    print(f"   exp_calls for ln({val_ln_3}) = ceil(log2({val_ln_3} / {PRECISION})) = {exp_calls_3}")
    print(f"   Total 'exp' calls (d): {exp_calls_1} + {exp_calls_2} + {exp_calls_3} = {result_d}")
    print("-" * 20)
    
    # --- Final Answer ---
    print("Final answer in a:b:c:d format:")
    print(f"{result_a}:{result_b}:{result_c}:{result_d}")


if __name__ == '__main__':
    solve_pioneer_problem()