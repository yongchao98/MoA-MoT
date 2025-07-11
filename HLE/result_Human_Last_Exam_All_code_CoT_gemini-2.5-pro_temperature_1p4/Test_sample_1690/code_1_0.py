import math

def solve_pioneer_problem():
    """
    Calculates the answers for the Pioneer probe problem based on the Wuxing architecture constraints.
    """
    # Given parameters
    P = 100  # Initial power in Watts
    x = 4    # Control subsystem power in Watts
    y = 5    # Camera subsystem power in Watts
    z = 1    # Sensor subsystem power in Watts
    half_life = 400 # Days

    # --- Part a: Time to turn off the camera subsystem ---
    # The camera turns off when the power is no longer sufficient for all three systems.
    # The decay equation is P_target = P * (1/2)^(t/half_life)
    # Solving for t gives: t = half_life * log2(P / P_target)
    P_target_camera = x + y + z
    t_camera_exact = half_life * math.log2(P / P_target_camera)
    t_camera_rounded = round(t_camera_exact)
    a = t_camera_rounded

    # --- Part b: Time to turn off the sensor subsystem ---
    # The sensor turns off when power is insufficient for the control and sensor systems.
    P_target_sensor = x + z
    t_sensor_exact = half_life * math.log2(P / P_target_sensor)
    t_sensor_rounded = round(t_sensor_exact)
    b = t_sensor_rounded

    # --- Part c: Memory usage in D for variables ---
    # A time- and memory-efficient program would reuse ln calculations.
    # The values needed are ln(100), ln(10), ln(5), ln(2). All can be derived from ln(2) and ln(5).
    # Minimal variables needed to be stored:
    # 4 ints for inputs (P, x, y, z): 4 * 5D = 20D
    # 2 fracs for calculated constants (ln(2), ln(5)): 2 * 6D = 12D
    # 2 ints for final integer results: 2 * 5D = 10D
    # Total memory usage = 20 + 12 + 10 = 42D
    mem_usage_inputs = 4 * 5
    mem_usage_ln_vals = 2 * 6
    mem_usage_results = 2 * 5
    c = mem_usage_inputs + mem_usage_ln_vals + mem_usage_results

    # --- Part d: Number of time program calls function exp ---
    # A time-efficient program only calls ln() for base values, ln(2) and ln(5).
    # The number of exp() calls in ln(v) is ceil(log2(v * 1e9)).
    exp_calls_for_ln2 = math.ceil(math.log2(2 * 1e9))
    exp_calls_for_ln5 = math.ceil(math.log2(5 * 1e9))
    d = exp_calls_for_ln2 + exp_calls_for_ln5

    # --- Print detailed explanations and final answer ---
    print(f"a. When to turn off the camera subsystem (days)?")
    print(f"   Equation: round({half_life} * log2({P} / ({x}+{y}+{z})))")
    print(f"   Calculation: round({half_life} * log2({P} / {P_target_camera})) = {a}\n")

    print(f"b. When to turn off the sensor subsystem (days)?")
    print(f"   Equation: round({half_life} * log2({P} / ({x}+{z})))")
    print(f"   Calculation: round({half_life} * log2({P} / {P_target_sensor})) = {b}\n")

    print(f"c. What is the memory usage in D for variables of this program?")
    print(f"   Calculation: {mem_usage_inputs}D (inputs) + {mem_usage_ln_vals}D (ln consts) + {mem_usage_results}D (results) = {c}D\n")
    
    print(f"d. What is the number of time this program call function exp?")
    print(f"   Calculation: calls for ln(2) + calls for ln(5) = {exp_calls_for_ln2} + {exp_calls_for_ln5} = {d}\n")
    
    # Final answer in the required format
    final_answer = f"{a}:{b}:{c}:{d}"
    print(f"<<<{final_answer}>>>")

solve_pioneer_problem()