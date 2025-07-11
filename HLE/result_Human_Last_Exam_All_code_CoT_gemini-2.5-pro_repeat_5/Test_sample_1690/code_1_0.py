import math

def solve_pioneer_problem():
    """
    Calculates the shutdown times and resource usage for the Pioneer probe
    based on the Wuxing architecture specifications.
    """
    # Input values from the problem
    P = 100  # Initial power in Watts
    x = 4    # Control subsystem power in Watts
    y = 5    # Camera subsystem power in Watts
    z = 1    # Sensor subsystem power in Watts
    T_half = 400 # Power half-life in days

    # a. Calculate when to turn off the camera subsystem
    # This happens when power is no longer sufficient for all three systems.
    # We solve for t in: P * (1/2)^(t/T_half) = x + y + z
    # t = T_half * log2(P / (x + y + z))
    power_req_camera = x + y + z
    time_camera_off = T_half * math.log2(P / power_req_camera)
    # The result should be rounded to the nearest total day
    answer_a = round(time_camera_off)

    # b. Calculate when to turn off the sensor subsystem
    # This happens when power is no longer sufficient for control and sensor systems.
    # We solve for t in: P * (1/2)^(t/T_half) = x + z
    # t = T_half * log2(P / (x + z))
    power_req_sensor = x + z
    time_sensor_off = T_half * math.log2(P / power_req_sensor)
    # The result should be rounded to the nearest total day
    answer_b = round(time_sensor_off)

    # c. Calculate the memory usage in D for the program's variables
    # An efficient Wuxing C program would need variables for:
    # P, x, y, z, T_half, ln_two, a reusable P_req, and a reusable result 't'.
    # This amounts to 8 'frac' variables.
    # Each 'frac' struct (signed char, unsigned char, signed char) uses
    # 2D + 2D + 2D = 6D of memory.
    num_frac_variables = 8
    memory_per_frac = 6 # in Decimal positions (D)
    answer_c = num_frac_variables * memory_per_frac

    # d. Calculate the number of times the program calls the 'exp' function
    # The formula t = T_half * (ln(P) - ln(P_req)) / ln(2) requires 4 ln calls:
    # ln(2), ln(P), ln(power_req_camera), and ln(power_req_sensor)
    # The number of exp calls for ln(v) is ceil(log2(v * 1e9))
    def count_exp_calls(value):
      return math.ceil(math.log2(value * 1e9))

    ln_call_values = [2, P, power_req_camera, power_req_sensor]
    answer_d = sum(count_exp_calls(v) for v in ln_call_values)

    # Print the final answer in the required a:b:c:d format
    print(f"{answer_a}:{answer_b}:{answer_c}:{answer_d}")

solve_pioneer_problem()
<<<1329:1729:48:135>>>