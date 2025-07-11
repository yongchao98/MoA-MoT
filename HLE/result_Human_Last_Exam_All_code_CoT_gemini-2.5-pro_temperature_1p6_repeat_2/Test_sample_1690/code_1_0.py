import math

def solve_pioneer_problem():
    """
    Calculates the answers to the Pioneer probe problem and its analysis on the Wuxing architecture.
    """
    # Input values
    P = 100
    x = 4
    y = 5
    z = 1

    # --- Part a: Time to turn off the camera subsystem ---
    # Power requirement with all subsystems on
    power_req_a = x + y + z
    # t = 400 * ln(P / P_req) / ln(2)
    time_a = 400 * math.log(P / power_req_a) / math.log(2)
    # The result should be rounded to total days
    answer_a = round(time_a)

    # --- Part b: Time to turn off the sensor subsystem ---
    # Power requirement with control and sensor subsystems on
    power_req_b = x + z
    time_b = 400 * math.log(P / power_req_b) / math.log(2)
    # The result should be rounded to total days
    answer_b = round(time_b)

    # --- Part c: Memory usage in D ---
    # An efficient C program would require variables for:
    # P, x, y, z (inputs)
    # ln2 (reused constant)
    # temp_res (reused for results)
    # Total = 6 variables
    num_frac_variables = 6
    # Each 'char' type is 2D. frac has 3 chars.
    memory_per_frac = 2 + 2 + 2  # D per frac
    answer_c = num_frac_variables * memory_per_frac

    # --- Part d: Number of 'exp' function calls ---
    # The program needs to calculate ln(2), ln(10), and ln(20).
    # The number of 'exp' calls for ln(x) is ceil(log2(x * 10^9)).
    
    # Value for part a's ln calculation
    val_a = P / power_req_a
    
    # Value for part b's ln calculation
    val_b = P / power_req_b
    
    calls_for_ln2 = math.ceil(math.log2(2 * 1e9))
    calls_for_lna = math.ceil(math.log2(val_a * 1e9))
    calls_for_lnb = math.ceil(math.log2(val_b * 1e9))
    
    answer_d = calls_for_ln2 + calls_for_lna + calls_for_lnb

    # Print the final result in the format a:b:c:d
    print(f"{answer_a}:{answer_b}:{answer_c}:{answer_d}")

solve_pioneer_problem()