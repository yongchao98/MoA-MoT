import math

def solve_pioneer_problem():
    """
    Solves the Pioneer probe power problem based on the Wuxing architecture specifications.
    """

    # --- Given Input Values ---
    P_initial = 100  # Initial power in Watts
    x = 4          # Control subsystem power in Watts
    y = 5          # Camera subsystem power in Watts
    z = 1          # Sensor subsystem power in Watts
    half_life_days = 400
    ln_precision = 1e-9

    # a. Calculate when the camera subsystem turns off
    power_req_cam = x + y + z
    # Equation for time t: P_req = P_initial * (0.5)^(t / half_life_days)
    # Solving for t: t = half_life_days * ln(P_req / P_initial) / ln(0.5)
    ratio_cam = power_req_cam / P_initial
    t_cam = half_life_days * (math.log(ratio_cam) / math.log(0.5))
    t_cam_rounded = int(round(t_cam))

    # b. Calculate when the sensor subsystem turns off
    power_req_sensor = x + z
    ratio_sensor = power_req_sensor / P_initial
    t_sensor = half_life_days * (math.log(ratio_sensor) / math.log(0.5))
    t_sensor_rounded = int(round(t_sensor))

    # c. Calculate the memory usage in D for an efficient C program
    # The most memory-efficient program would use the smallest types:
    # - frac P = 100;          (6D)
    # - frac four_hundred = 400; (6D)
    # - char x=4, y=5, z=1;   (2D * 3 = 6D)
    # - frac ln_half;            (6D, stores result of ln(0.5) for reuse)
    # - frac t_result;         (6D, reused for both calculations)
    # - int days_cam;           (5D, result is large)
    # - int days_sensor;        (5D, result is large)
    # Total memory usage = 6 + 6 + (3 * 2) + 6 + 6 + 5 + 5
    memory_usage = 40

    # d. Calculate the number of times the 'exp' function is called
    # The 'ln(x)' function calls 'exp' once per loop iteration.
    # The number of iterations 'k' is determined by `ceil(log2(x / precision))`.
    def get_exp_calls_for_ln(val, precision):
      # For ln(val), initial interval is [0, val].
      # Number of iterations is k = ceil(log2(val / precision))
      return math.ceil(math.log2(val / precision))

    # An efficient program makes 3 ln() calls: ln(0.5), ln(ratio_cam), ln(ratio_sensor)
    exp_calls_for_ln_half = get_exp_calls_for_ln(0.5, ln_precision)
    exp_calls_for_ln_cam = get_exp_calls_for_ln(ratio_cam, ln_precision)
    exp_calls_for_ln_sensor = get_exp_calls_for_ln(ratio_sensor, ln_precision)
    
    total_exp_calls = exp_calls_for_ln_half + exp_calls_for_ln_cam + exp_calls_for_ln_sensor

    # Print the final result in the format a:b:c:d
    print(f"{t_cam_rounded}:{t_sensor_rounded}:{memory_usage}:{total_exp_calls}")

solve_pioneer_problem()
<<<1329:1729:40:82>>>