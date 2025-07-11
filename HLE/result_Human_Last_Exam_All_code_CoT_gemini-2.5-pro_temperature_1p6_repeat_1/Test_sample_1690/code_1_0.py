import math

def solve_pioneer_power():
    """
    Calculates the answers for the Pioneer probe power management problem.
    """
    # Given parameters
    P = 100.0
    x = 4.0
    y = 5.0
    z = 1.0
    T_half = 400.0

    # --- Part a: Camera turn-off time ---
    power_target_cam = x + y + z
    # The formula is: t = T_half * ln(P / power_target) / ln(2)
    t_cam = T_half * math.log(P / power_target_cam) / math.log(2)
    a = round(t_cam)

    print(f"Equation for camera turn-off (a): t = {T_half} * ln({P} / {power_target_cam}) / ln(2)")
    
    # --- Part b: Sensor turn-off time ---
    power_target_sensor = x + z
    # The formula is the same, just with a new target power
    t_sensor = T_half * math.log(P / power_target_sensor) / math.log(2)
    b = round(t_sensor)

    print(f"Equation for sensor turn-off (b): t = {T_half} * ln({P} / {power_target_sensor}) / ln(2)")

    # --- Part c: Memory usage in D ---
    # A memory-efficient C program would reuse variables. The minimal set of `frac`
    # variables is for: T_half, P, target_power, two, ln2, result.
    # Each frac is 6D (2D for n, 2D for d, 2D for e).
    num_frac_variables = 6
    memory_per_frac = 6
    c = num_frac_variables * memory_per_frac

    # --- Part d: Number of exp() calls ---
    # A time-efficient C program calculates ln(2) and ln(10) once, then computes
    # ln(20) as ln(2) + ln(10). The ln(x) implementation calls exp() in a loop.
    # The number of calls is approximately ceil(log2(x) + 9*log2(10)).
    # log2(10) is approx 3.3219. The constant part is 9*log2(10) approx 29.897.
    exp_calls_for_ln2 = math.ceil(math.log2(2) + 9 * math.log2(10)) # 1 + 29.897 -> 31
    exp_calls_for_ln10 = math.ceil(math.log2(10) + 9 * math.log2(10)) # 3.32 + 29.897 -> 34
    d = exp_calls_for_ln2 + exp_calls_for_ln10
    
    # --- Final Answer ---
    final_answer = f"{a}:{b}:{c}:{d}"
    print(f"\nFinal calculated values:")
    print(f"a = {a} days")
    print(f"b = {b} days")
    print(f"c = {c} D")
    print(f"d = {d} calls")

    print(f"\n<<<{final_answer}>>>")

solve_pioneer_power()