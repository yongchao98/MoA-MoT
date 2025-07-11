import math

def solve_pioneer_problem():
    """
    Calculates the answers for the Pioneer probe problem and explains the steps.
    """
    # Given values
    P = 100.0
    x = 4.0
    y = 5.0
    z = 1.0

    # --- Part a & b: Calculate Shutdown Times ---
    
    # The governing equation for power decay is P(t) = P * (1/2)^(t/400).
    # Solving for t gives the formula: t = 400 * ln(P / P_req) / ln(2).

    # a. Camera subsystem shutdown time
    P_req_cam = x + y + z
    t_cam_off = 400 * math.log(P / P_req_cam) / math.log(2)
    t_cam_off_rounded = int(round(t_cam_off))

    # b. Sensor subsystem shutdown time
    P_req_sensor = x + z
    t_sensor_off = 400 * math.log(P / P_req_sensor) / math.log(2)
    t_sensor_off_rounded = int(round(t_sensor_off))

    # --- Part c: Calculate Memory Usage ---

    # To be time-efficient, the program should minimize calls to ln(), the most
    # expensive function besides exp().
    # t_cam_off needs ln(100/10) = ln(10)
    # t_sensor_off needs ln(100/5) = ln(20)
    # The optimization is ln(20) = ln(2*10) = ln(2) + ln(10).
    # This requires calling ln() only twice (for 2 and 10).
    #
    # Variables needed in main():
    # 4 for inputs: P, x, y, z
    # 2 for intermediate ln values: ln_2, ln_10
    # 2 for results: t_cam_off, t_sensor_off
    # Total variables = 8.
    num_frac_variables = 8
    # Memory for one 'frac' variable = sizeof(n) + sizeof(d) + sizeof(e)
    # = 2D + 2D + 2D = 6D.
    mem_usage = num_frac_variables * 6

    # --- Part d: Calculate Number of exp() Calls ---

    # The number of iterations 'k' in ln(x) is ceil(log2(x * 10^9)).
    # Each iteration calls exp() once.
    # The efficient program calls ln(2) and ln(10).
    
    log2_1e9 = math.log2(1e9)
    # Calls for ln(2)
    exp_calls_ln2 = math.ceil(math.log2(2) + log2_1e9)
    # Calls for ln(10)
    exp_calls_ln10 = math.ceil(math.log2(10) + log2_1e9)
    total_exp_calls = exp_calls_ln2 + exp_calls_ln10

    # --- Print Detailed Explanation and Results ---

    print("Here is the step-by-step calculation:")

    print("\na. When to turn off the camera subsystem (days)?")
    print(f"The camera is turned off when power drops below the requirement for all systems: x+y+z = {x:.0f}+{y:.0f}+{z:.0f} = {P_req_cam:.0f} W.")
    print("The final equation for time 't' is: t = 400 * ln(P / P_req) / ln(2)")
    print(f"t = 400 * ln({P:.0f} / {P_req_cam:.0f}) / ln(2)")
    print(f"t = 400 * ln({P/P_req_cam:.0f}) / ln(2)")
    print(f"t = 400 * {math.log(P/P_req_cam):.6f} / {math.log(2):.6f} = {t_cam_off:.2f}")
    print(f"Rounded to the nearest day, the answer is: {t_cam_off_rounded}")

    print("\nb. When to turn off the sensor subsystem (days)?")
    print(f"The sensor is turned off when power drops below the requirement for the control and sensor systems: x+z = {x:.0f}+{z:.0f} = {P_req_sensor:.0f} W.")
    print(f"t = 400 * ln({P:.0f} / {P_req_sensor:.0f}) / ln(2)")
    print(f"t = 400 * ln({P/P_req_sensor:.0f}) / ln(2)")
    print(f"t = 400 * {math.log(P/P_req_sensor):.6f} / {math.log(2):.6f} = {t_sensor_off:.2f}")
    print(f"Rounded to the nearest day, the answer is: {t_sensor_off_rounded}")
    
    print("\nc. What is the memory usage in D for variables of this program?")
    print(f"A time-efficient program uses {num_frac_variables} 'frac' variables in main().")
    print(f"Each 'frac' is 6D (3 chars * 2D/char).")
    print(f"Total memory = {num_frac_variables} * 6D = {mem_usage}D.")

    print("\nd. What is the number of time this program call function exp?")
    print("A time-efficient program calls ln(2) and ln(10). Each ln() call iterates and calls exp().")
    print(f"Number of exp() calls for ln(2) = ceil(log2(2*10^9)) = {exp_calls_ln2}")
    print(f"Number of exp() calls for ln(10) = ceil(log2(10*10^9)) = {exp_calls_ln10}")
    print(f"Total exp() calls = {exp_calls_ln2} + {exp_calls_ln10} = {total_exp_calls}")

    # --- Final Answer ---
    final_answer = f"{t_cam_off_rounded}:{t_sensor_off_rounded}:{mem_usage}:{total_exp_calls}"
    print("\nFinal Answer (a:b:c:d):")
    print(f"<<<{final_answer}>>>")

solve_pioneer_problem()