import math

def solve_wuxing_problem():
    """
    Calculates the solutions for the Pioneer probe problem based on the Wuxing architecture.
    """
    # --- Problem Parameters ---
    P_initial = 100.0
    x_power = 4.0
    y_power = 5.0
    z_power = 1.0
    T_half = 400.0

    # --- Part a: Camera subsystem shutdown time ---
    # Shutdown occurs when power is no longer sufficient for all subsystems.
    # Final Power (P_final) = x + y + z
    # Equation: t = T_half * (ln(P_initial) - ln(P_final)) / ln(2)
    power_req_cam_off = x_power + y_power + z_power
    t_cam_exact = T_half * (math.log(P_initial) - math.log(power_req_cam_off)) / math.log(2)
    a_t_cam_rounded = int(round(t_cam_exact))

    # --- Part b: Sensor subsystem shutdown time ---
    # Shutdown occurs when power is no longer sufficient for control and sensor.
    # Final Power (P_final) = x + z
    # Equation: t = T_half * (ln(P_initial) - ln(P_final)) / ln(2)
    power_req_sensor_off = x_power + z_power
    t_sensor_exact = T_half * (math.log(P_initial) - math.log(power_req_sensor_off)) / math.log(2)
    b_t_sensor_rounded = int(round(t_sensor_exact))

    # --- Part c: Memory usage in D ---
    # An efficient C program would store inputs, calculation variables, and results.
    # It would cache ln(P) to be time-efficient.
    # int size = 5D, frac size = 6D
    # Variables: 4 ints (P,x,y,z), 6 fracs (P_f, p_final, t_days, T_half, ln2, lnP), 2 ints (results)
    # Final Equation: 4*5 + 6*6 + 2*5
    num_input_ints = 4
    size_int = 5
    num_frac_vars = 6
    size_frac = 6
    num_result_ints = 2
    c_memory_usage = (num_input_ints * size_int) + (num_frac_vars * size_frac) + (num_result_ints * size_int)

    # --- Part d: Number of calls to exp() ---
    # The ln(x) function calls exp() roughly ceil(log2(x) + 30) times.
    # A time-efficient program computes ln(2), ln(100), ln(10), and ln(5) once each.
    def count_exp_calls(val):
        # Number of iterations is ceil(log2(initial_range / precision))
        # ceil(log2(val / 1e-9)) = ceil(log2(val) + 9*log2(10)) ~ ceil(log2(val)+30)
        if val <= 0: return 0
        return math.ceil(math.log2(val) + 30)
    
    # Values for which ln() is computed in an efficient program:
    # ln(2)
    # ln(P_initial) which is ln(100)
    # ln(power_req_cam_off) which is ln(10)
    # ln(power_req_sensor_off) which is ln(5)
    
    calls_for_ln2 = count_exp_calls(2.0)
    calls_for_ln_P = count_exp_calls(P_initial)
    calls_for_ln_cam = count_exp_calls(power_req_cam_off)
    calls_for_ln_sensor = count_exp_calls(power_req_sensor_off)

    d_total_exp_calls = calls_for_ln2 + calls_for_ln_P + calls_for_ln_cam + calls_for_ln_sensor

    # --- Final Output ---
    # The prompt asks to output each number in the final equation. 
    # The calculations above show how each final number is derived.
    # For a: t = 400 * (ln(100) - ln(10)) / ln(2)
    # For b: t = 400 * (ln(100) - ln(5)) / ln(2)
    # For c: mem = 4*5 + 6*6 + 2*5
    # For d: calls = calls(2) + calls(100) + calls(10) + calls(5)
    print(f"{a_t_cam_rounded}:{b_t_sensor_rounded}:{c_memory_usage}:{d_total_exp_calls}")

solve_wuxing_problem()