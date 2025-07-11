import math

def solve_pioneer_problem():
    """
    Calculates the answers to the Pioneer probe problem based on the Wuxing architecture specifications.
    """
    # --- Input values from the problem ---
    P = 100.0  # Initial power
    x = 4.0    # Control subsystem power
    y = 5.0    # Camera subsystem power
    z = 1.0    # Sensor subsystem power
    half_life = 400.0 # days
    precision = 1e-9

    # --- Calculations ---

    # a. When the camera subsystem is turned off
    power_total = x + y + z
    # t = 400 * (ln(100) - ln(10)) / ln(2) = 400 * ln(10) / ln(2)
    t_camera_exact = half_life * (math.log(P) - math.log(power_total)) / math.log(2)
    a_camera_shutdown_days = round(t_camera_exact)

    # b. When the sensor subsystem is turned off
    power_sensor_only = x + z
    # t = 400 * (ln(100) - ln(5)) / ln(2) = 400 * ln(20) / ln(2)
    t_sensor_exact = half_life * (math.log(P) - math.log(power_sensor_only)) / math.log(2)
    b_sensor_shutdown_days = round(t_sensor_exact)

    # c. Memory usage in D for variables
    # An efficient C program would need variables for:
    # - Inputs: P, x, y, z (4 frac)
    # - A reusable constant for 1/ln(2) (1 frac)
    # - Two output variables: t_camera, t_sensor (2 frac)
    # Total variables = 4 + 1 + 2 = 7
    # Size of one frac = signed char (2D) + unsigned char (2D) + signed char (2D) = 6D
    num_variables = 7
    size_of_frac_in_D = 6
    c_memory_usage_D = num_variables * size_of_frac_in_D

    # d. Number of times the program calls exp()
    # An efficient program calls ln() 4 times: ln(2), ln(5), ln(10), ln(100)
    # The number of exp() calls in ln(v) is ceil(log2(v / precision))
    ln_inputs = [2.0, 5.0, 10.0, 100.0]
    total_exp_calls = 0
    for val in ln_inputs:
        # Number of iterations k satisfies: v / (2^k) <= precision => 2^k >= v / precision
        # k >= log2(v / precision)
        k = math.ceil(math.log2(val / precision))
        total_exp_calls += k
    d_exp_calls = total_exp_calls
    
    # --- Final Output ---
    # The prompt asks to output each number in the final equation.
    # Here are the numbers that form the final answer string.
    print(f"a (Camera shutdown day): {a_camera_shutdown_days}")
    print(f"b (Sensor shutdown day): {b_sensor_shutdown_days}")
    print(f"c (Memory usage in D): {c_memory_usage_D}")
    print(f"d (Number of exp calls): {d_exp_calls}")
    
    # Combine into the final answer format
    final_answer = f"{int(a_camera_shutdown_days)}:{int(b_sensor_shutdown_days)}:{int(c_memory_usage_D)}:{int(d_exp_calls)}"
    print("\nFinal Answer String:")
    print(f"<<<{final_answer}>>>")


solve_pioneer_problem()