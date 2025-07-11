import math

def solve_wuxing_problem():
    """
    Calculates the answers for the Pioneer probe problem based on the Wuxing architecture.
    """
    
    # --- Part a & b: Calculate shutdown times ---
    P = 100.0
    x = 4.0
    y = 5.0
    z = 1.0
    T_half = 400.0

    # a. Time to turn off the camera subsystem
    # This happens when power is no longer sufficient for control + camera + sensor.
    # The equation is t_a = T_half * log2(P / (x + y + z))
    required_power_a = x + y + z
    t_a = T_half * math.log2(P / required_power_a)
    answer_a = int(round(t_a))

    # b. Time to turn off the sensor subsystem
    # This happens when power is no longer sufficient for control + sensor.
    # The equation is t_b = T_half * log2(P / (x + z))
    required_power_b = x + z
    t_b = T_half * math.log2(P / required_power_b)
    answer_b = int(round(t_b))

    # --- Part c: Calculate memory usage ---
    # Memory for variables in the most efficient Wuxing C program.
    # P: int (5D) as it is > 99
    # x, y, z: digit (1D each)
    # 3 frac variables for calculations (result, ln_P, ln_two), each 6D.
    # Total D = 5 + (3 * 1) + (3 * 6)
    mem_p = 5
    mem_xyz = 3 * 1
    mem_frac_vars = 3 * 6
    answer_c = mem_p + mem_xyz + mem_frac_vars

    # --- Part d: Calculate the number of exp function calls ---
    # The program needs to calculate ln(2), ln(100), ln(10), and ln(5).
    # The number of exp calls for ln(x) is k = ceil(log2(x) + 9*log2(10)).
    
    def count_exp_calls(val):
        """Calculates the number of iterations (exp calls) for ln(val)."""
        # A small epsilon is subtracted to handle potential floating point inaccuracies before ceiling.
        log2_term = math.log2(val)
        log2_10_pow_9 = 9 * math.log2(10)
        return math.ceil(log2_term + log2_10_pow_9 - 1e-9)

    calls_for_ln_2 = count_exp_calls(2)    # ln(2)
    calls_for_ln_10 = count_exp_calls(10)   # ln(10) = ln(x+y+z)
    calls_for_ln_5 = count_exp_calls(5)     # ln(5) = ln(x+z)
    calls_for_ln_100 = count_exp_calls(100)  # ln(100) = ln(P)
    
    answer_d = calls_for_ln_2 + calls_for_ln_5 + calls_for_ln_10 + calls_for_ln_100

    # Print the final answer in the format a:b:c:d
    print(f"{answer_a}:{answer_b}:{answer_c}:{answer_d}")

solve_wuxing_problem()
<<<1329:1729:26:135>>>