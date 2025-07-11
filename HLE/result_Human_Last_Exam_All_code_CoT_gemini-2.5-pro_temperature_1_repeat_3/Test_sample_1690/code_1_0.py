import math

def solve_pioneer_problem():
    """
    Calculates the answers for the Pioneer probe problem based on the Wuxing architecture.
    """
    # --- Problem Setup ---
    P = 100  # Initial power (W)
    x = 4    # Control subsystem power (W)
    y = 5    # Camera subsystem power (W)
    z = 1    # Sensor subsystem power (W)
    T_half = 400 # Power half-life (days)
    
    print("--- Solving the Pioneer Probe Power Problem ---")
    print(f"Initial Values: P={P}W, x={x}W, y={y}W, z={z}W, T_half={T_half} days\n")

    # --- Part a: Camera Shutdown Time ---
    print("a. When to turn off the camera subsystem?")
    # This occurs when power P(t) drops below the total required power (x + y + z).
    # The equation is: P * (1/2)^(t / T_half) = (x + y + z)
    # Solving for t: t = T_half * ln(P / (x + y + z)) / ln(2)
    req_power_all = x + y + z
    ratio_all = P / req_power_all
    t_cam_off_exact = T_half * math.log(ratio_all) / math.log(2)
    a = round(t_cam_off_exact)

    print(f"Required power for all systems (x+y+z) = {x} + {y} + {z} = {req_power_all} W")
    print(f"Equation: t = {T_half} * ln({P} / {req_power_all}) / ln(2)")
    print(f"Calculation: t = {T_half} * ln({ratio_all:.2f}) / ln(2) = {t_cam_off_exact:.2f} days")
    print(f"Rounded to the nearest day: {a} days\n")

    # --- Part b: Sensor Shutdown Time ---
    print("b. When to turn off the sensor subsystem?")
    # This occurs when power P(t) drops below the remaining required power (x + z).
    # The equation is: P * (1/2)^(t / T_half) = (x + z)
    # Solving for t: t = T_half * ln(P / (x + z)) / ln(2)
    req_power_no_cam = x + z
    ratio_no_cam = P / req_power_no_cam
    t_sensor_off_exact = T_half * math.log(ratio_no_cam) / math.log(2)
    b = round(t_sensor_off_exact)

    print(f"Required power without camera (x+z) = {x} + {z} = {req_power_no_cam} W")
    print(f"Equation: t = {T_half} * ln({P} / {req_power_no_cam}) / ln(2)")
    print(f"Calculation: t = {T_half} * ln({ratio_no_cam:.2f}) / ln(2) = {t_sensor_off_exact:.2f} days")
    print(f"Rounded to the nearest day: {b} days\n")

    # --- Part c: Memory Usage ---
    print("c. What is the memory usage in D for the most efficient program?")
    # Wuxing data type sizes: char=2D, int=5D, frac=6D.
    # Minimal variable set for a memory-efficient C program:
    #   - 4 'char' variables for inputs P,x,y,z.
    #   - 3 'frac' variables (for ln(2), the ratio, and the time calculation), reused.
    #   - 2 'int' variables for the final rounded day counts.
    size_char = 2
    size_int = 5
    size_frac = 6
    num_chars = 4
    num_fracs = 3
    num_ints = 2
    c = (num_chars * size_char) + (num_fracs * size_frac) + (num_ints * size_int)
    print(f"Variables: {num_chars}*char ({size_char}D) + {num_fracs}*frac ({size_frac}D) + {num_ints}*int ({size_int}D)")
    print(f"Calculation: ({num_chars} * {size_char}) + ({num_fracs} * {size_frac}) + ({num_ints} * {size_int}) = {c} D\n")

    # --- Part d: Number of 'exp' function calls ---
    print("d. How many times is the 'exp' function called?")
    # The 'exp' function is called once per iteration of the 'ln' function's loop.
    # Iterations for ln(val) is ceil(log2(val / epsilon)), with epsilon=1e-9.
    # The most time-efficient program makes 3 calls to ln(): ln(10), ln(20), and ln(2).
    def ln_iterations(val):
        epsilon = 1e-9
        return math.ceil(math.log2(val / epsilon))

    calls_for_ln_ratio_all = ln_iterations(ratio_all)
    calls_for_ln_ratio_no_cam = ln_iterations(ratio_no_cam)
    calls_for_ln_2 = ln_iterations(2)
    d = calls_for_ln_ratio_all + calls_for_ln_ratio_no_cam + calls_for_ln_2
    
    print(f"Total calls = calls for ln({ratio_all}) + calls for ln({ratio_no_cam}) + calls for ln(2)")
    print(f"Calculation: {calls_for_ln_ratio_all} + {calls_for_ln_ratio_no_cam} + {calls_for_ln_2} = {d}\n")
    
    # --- Final Answer ---
    print("--- Final Answer ---")
    print("The final answer in the format a:b:c:d is:")
    print(f"<<<{a}:{b}:{c}:{d}>>>")

solve_pioneer_problem()