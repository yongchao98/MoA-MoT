import math

def solve_pioneer_problem():
    """
    Solves the Pioneer probe power problem based on the Wuxing architecture specifications.
    """
    # Given values
    P = 100.0  # Initial power in Watts
    x = 4.0    # Control subsystem power in Watts
    y = 5.0    # Camera subsystem power in Watts
    z = 1.0    # Sensor subsystem power in Watts
    T_half = 400.0  # Half-life in days

    # --- Part a & b: Calculate shutdown times ---

    # Power level when camera is turned off
    power_target_camera = x + y + z
    # Power level when sensor is turned off
    power_target_sensor = x + z

    # The formula to find the time 't' is derived from P(t) = P_initial * (1/2)^(t/T_half)
    # t = -T_half * ln(P_target / P_initial) / ln(2)

    # a. Time to turn off the camera subsystem
    t_camera = -T_half * math.log(power_target_camera / P) / math.log(2)
    days_camera_off = int(round(t_camera))

    # b. Time to turn off the sensor subsystem
    t_sensor = -T_half * math.log(power_target_sensor / P) / math.log(2)
    days_sensor_off = int(round(t_sensor))

    # --- Part c: Calculate Memory Usage ---

    # For a memory-efficient program, we need variables for:
    # Inputs: P, x, y, z
    # Constant: T_half
    # Intermediates for efficiency: ln(P), ln(2)
    # Outputs: t_cam, t_sensor
    # Total variables = 4 (inputs) + 1 (const) + 2 (intermediate) + 2 (output) = 9
    # Note: In a real C program, inputs might be read directly, but they still occupy
    # memory when loaded into variables for calculation. We count all declared variables in main.
    num_frac_variables = 9
    
    # Each 'frac' has 3 'char' members (n, d, e)
    chars_per_frac = 3
    # Each 'char' on Wuxing is 2 decimal digits (D)
    d_per_char = 2
    
    # c. Total memory usage in D
    memory_usage_d = num_frac_variables * chars_per_frac * d_per_char

    # --- Part d: Calculate number of 'exp' function calls ---

    # The number of calls to exp() inside ln(x) is ceil(log2(x / 1e-9))
    def count_exp_calls(val):
        if val <= 0:
            return 0
        # The while loop condition is (r - l > 1e-9), where the initial interval (r-l) is x.
        # The number of binary search steps k must satisfy x / (2^k) <= 1e-9
        # k >= log2(x / 1e-9)
        precision = 1e-9
        return math.ceil(math.log2(val / precision))

    # The efficient program calculates ln for 4 unique values: P, 2, power_target_camera, power_target_sensor
    calls_for_ln_P = count_exp_calls(P)
    calls_for_ln_2 = count_exp_calls(2)
    calls_for_ln_cam = count_exp_calls(power_target_camera)
    calls_for_ln_sensor = count_exp_calls(power_target_sensor)
    
    # d. Total number of exp calls
    total_exp_calls = calls_for_ln_P + calls_for_ln_2 + calls_for_ln_cam + calls_for_ln_sensor

    # --- Final Output ---
    # Print the results in the format a:b:c:d
    print(f"{days_camera_off}:{days_sensor_off}:{memory_usage_d}:{total_exp_calls}")

solve_pioneer_problem()