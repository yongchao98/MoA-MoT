import math

def solve_wuxing_problem():
    """
    Calculates the answers for the Wuxing probe problem and prints the steps.
    """
    # Given initial values
    P = 100.0
    x = 4.0
    y = 5.0
    z = 1.0
    T_half = 400.0

    # --- Part a & b: Shutdown Times ---
    print("--- Calculating Shutdown Times (a, b) ---")
    
    # Power requirements
    P_req_a = x + y + z  # For all systems
    P_req_b = x + z      # For control and sensor only

    # Using the formula: t = T_half * (ln(P) - ln(P_req)) / ln(2)
    # Python's math.log is the natural logarithm (ln)
    
    # Calculate values for the equation
    ln_P = math.log(P)
    ln_req_a = math.log(P_req_a)
    ln_req_b = math.log(P_req_b)
    ln_2 = math.log(2)

    # a. Camera subsystem shutdown time
    t_a_float = T_half * (ln_P - ln_req_a) / ln_2
    t_a_days = round(t_a_float)

    print("\na. Time to turn off the camera subsystem (when power < {} W):".format(P_req_a))
    print("t_a = T_half * (ln(P) - ln(x+y+z)) / ln(2)")
    print("{:.0f} days = {:.0f} * (ln({:.0f}) - ln({:.0f})) / ln(2)".format(t_a_days, T_half, P, P_req_a))
    print("{:.0f} days = {:.0f} * ({:.4f} - {:.4f}) / {:.4f}".format(t_a_days, T_half, ln_P, ln_req_a, ln_2))
    
    # b. Sensor subsystem shutdown time
    t_b_float = T_half * (ln_P - ln_req_b) / ln_2
    t_b_days = round(t_b_float)

    print("\nb. Time to turn off the sensor subsystem (when power < {} W):".format(P_req_b))
    print("t_b = T_half * (ln(P) - ln(x+z)) / ln(2)")
    print("{:.0f} days = {:.0f} * (ln({:.0f}) - ln({:.0f})) / ln(2)".format(t_b_days, T_half, P, P_req_b))
    print("{:.0f} days = {:.0f} * ({:.4f} - {:.4f}) / {:.4f}".format(t_b_days, T_half, ln_P, ln_req_b, ln_2))
    
    # --- Part c: Memory Usage ---
    print("\n--- Calculating Memory Usage (c) ---")
    
    # Size of one frac variable = signed char (2D) + unsigned char (2D) + signed char (2D) = 6D
    size_of_frac = 6 # in Decimal Digits (D)
    
    # Most efficient program needs to store results of key calculations.
    # We need variables for: ln(2), ln(100), result t_a, and result t_b.
    num_variables = 4
    memory_usage = num_variables * size_of_frac
    
    print("Memory for program variables in an efficient C implementation:")
    print("The program must store at least 4 frac values: ln(2), ln(100), t_a, and t_b.")
    print("Size of frac = 2D (n) + 2D (d) + 2D (e) = 6D")
    print("Total memory usage = {} variables * {} D/variable = {}D".format(num_variables, size_of_frac, memory_usage))

    # --- Part d: Number of exp() Calls ---
    print("\n--- Calculating exp() Function Calls (d) ---")
    
    # The program needs to compute ln(2), ln(10), ln(5), and ln(100).
    # Number of exp calls inside ln(val) is ceil(log2(val * 1e9))
    
    calls_ln2   = math.ceil(math.log2(2 * 1e9))
    calls_ln5   = math.ceil(math.log2(5 * 1e9))
    calls_ln10  = math.ceil(math.log2(10 * 1e9))
    calls_ln100 = math.ceil(math.log2(100 * 1e9))
    total_exp_calls = calls_ln2 + calls_ln5 + calls_ln10 + calls_ln100
    
    print("Total calls to exp() is the sum of calls from each ln() calculation.")
    print("Calls from ln(2):   ceil(log2(2 * 1e9))   = {}".format(calls_ln2))
    print("Calls from ln(5):   ceil(log2(5 * 1e9))   = {}".format(calls_ln5))
    print("Calls from ln(10):  ceil(log2(10 * 1e9))  = {}".format(calls_ln10))
    print("Calls from ln(100): ceil(log2(100 * 1e9)) = {}".format(calls_ln100))
    print("Total exp() calls = {} + {} + {} + {} = {}".format(calls_ln2, calls_ln5, calls_ln10, calls_ln100, total_exp_calls))

    # --- Final Answer ---
    print("\n--- Final Answer ---")
    final_answer = "{}:{}:{}:{}".format(t_a_days, t_b_days, memory_usage, total_exp_calls)
    print(final_answer)
    return final_answer

solve_wuxing_problem()