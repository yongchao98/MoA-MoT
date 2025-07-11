import math

def solve_pioneer_problem():
    """
    Calculates the shutdown times and resource usage for the Pioneer probe
    based on the Wuxing architecture constraints.
    """
    # Input values from the problem description
    P = 100.0
    x = 4.0
    y = 5.0
    z = 1.0
    T_half = 400.0

    # a. Calculate when the camera subsystem is turned off
    power_threshold_a = x + y + z
    # The formula is t = T_half * ln(P_initial / P_threshold) / ln(2)
    t_camera_float = T_half * math.log(P / power_threshold_a) / math.log(2)
    a = int(round(t_camera_float))

    # b. Calculate when the sensor subsystem is turned off
    power_threshold_b = x + z
    t_sensor_float = T_half * math.log(P / power_threshold_b) / math.log(2)
    b = int(round(t_sensor_float))

    # c. Calculate the memory usage in D for variables of the most efficient program
    # The program needs to store constants (400, 2, 10) and results (ln2, ln10, t_camera, t_sensor).
    # This is 7 `frac` variables in total.
    # Each frac is composed of 3 char types (2D each), so 3 * 2D = 6D.
    num_frac_variables = 7
    size_of_frac_in_D = 6
    c = num_frac_variables * size_of_frac_in_D

    # d. Calculate the number of times the expensive 'exp' function is called.
    # The most efficient program calls ln(2) and ln(10).
    # Number of exp calls for ln(v) is ceil(log2(v * 1e9)).
    # We use math.log2 for log base 2.
    precision_factor = 1e9
    exp_calls_for_ln2 = math.ceil(math.log2(2 * precision_factor))
    exp_calls_for_ln10 = math.ceil(math.log2(10 * precision_factor))
    d = exp_calls_for_ln2 + exp_calls_for_ln10

    # Print the final answer in the format a:b:c:d
    print(f"{a}:{b}:{c}:{d}")

solve_pioneer_problem()