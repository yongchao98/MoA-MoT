import math

def solve():
    """
    Solves the Pioneer probe problem based on the Wuxing architecture specifications.
    """
    # --- Input values ---
    P = 100.0  # Initial power
    x = 4.0   # Control subsystem power
    y = 5.0   # Camera subsystem power
    z = 1.0   # Sensor subsystem power
    H = 400.0  # Half-life in days

    print("Thinking Process & Calculations:")
    print("---------------------------------")
    print(f"Inputs: P={P}, x={x}, y={y}, z={z}, H={H}")
    print("The core formula to find the time 't' until power drops to a required level 'R' is:")
    print("t = H * (ln(P) - ln(R)) / ln(2)\n")

    # --- Part a: Camera shutdown time ---
    print("a. Calculating camera shutdown time:")
    R_cam = x + y + z
    t_cam = H * (math.log(P) - math.log(R_cam)) / math.log(2)
    t_cam_rounded = int(round(t_cam))
    print(f"   - Required power for all systems (R_cam) = {x} + {y} + {z} = {R_cam} W")
    print(f"   - t_cam = {H} * (ln({P}) - ln({R_cam})) / ln(2) = {t_cam:.2f} days")
    print(f"   - Rounded to the nearest day, a = {t_cam_rounded}\n")


    # --- Part b: Sensor shutdown time ---
    print("b. Calculating sensor shutdown time:")
    R_sensor = x + z
    t_sensor = H * (math.log(P) - math.log(R_sensor)) / math.log(2)
    t_sensor_rounded = int(round(t_sensor))
    print(f"   - Required power for control+sensor (R_sensor) = {x} + {z} = {R_sensor} W")
    print(f"   - t_sensor = {H} * (ln({P}) - ln({R_sensor})) / ln(2) = {t_sensor:.2f} days")
    print(f"   - Rounded to the nearest day, b = {t_sensor_rounded}\n")


    # --- Part c: Memory usage ---
    print("c. Calculating memory usage:")
    print("   - An efficient C program would declare variables for inputs, constants, and results.")
    print("   - Variables needed: P, x, y, z, H (5 inputs), common_factor (H/ln2), lnP, t_cam, t_sensor (4 for calculation/results).")
    num_frac_vars = 9
    size_of_frac_in_D = 6 # 3 chars * 2D/char
    mem_usage = num_frac_vars * size_of_frac_in_D
    print(f"   - Total variables: {num_frac_vars} of type 'frac'.")
    print(f"   - Size of one 'frac' variable = 3 * sizeof(char) = 3 * 2D = {size_of_frac_in_D}D.")
    print(f"   - Total memory usage = {num_frac_vars} * {size_of_frac_in_D}D = {mem_usage}D.")
    print(f"   - Therefore, c = {mem_usage}\n")


    # --- Part d: Number of exp() calls ---
    print("d. Calculating the number of exp() function calls:")
    print("   - The program calls ln() 4 times: ln(2), ln(P), ln(R_cam), ln(R_sensor).")
    print("   - Each ln(v) call makes k = ceil(log2(v * 10^9)) calls to exp().")

    def count_exp_calls(val):
        # k = ceil(log2(v * 10^9))
        return math.ceil(math.log2(val * 1e9))

    exp_calls_for_ln2 = count_exp_calls(2)
    exp_calls_for_lnP = count_exp_calls(P)
    exp_calls_for_lnR_cam = count_exp_calls(R_cam)
    exp_calls_for_lnR_sensor = count_exp_calls(R_sensor)
    total_exp_calls = exp_calls_for_ln2 + exp_calls_for_lnP + exp_calls_for_lnR_cam + exp_calls_for_lnR_sensor

    print(f"   - Calls for ln(2): {exp_calls_for_ln2}")
    print(f"   - Calls for ln({P}): {exp_calls_for_lnP}")
    print(f"   - Calls for ln({R_cam}): {exp_calls_for_lnR_cam}")
    print(f"   - Calls for ln({R_sensor}): {exp_calls_for_lnR_sensor}")
    print(f"   - Total exp() calls = {exp_calls_for_ln2} + {exp_calls_for_lnP} + {exp_calls_for_lnR_cam} + {exp_calls_for_lnR_sensor} = {total_exp_calls}")
    print(f"   - Therefore, d = {total_exp_calls}\n")
    
    # --- Final Answer ---
    final_answer = f"{t_cam_rounded}:{t_sensor_rounded}:{mem_usage}:{total_exp_calls}"
    print("---------------------------------")
    print(f"Final combined answer (a:b:c:d):")
    print(f"<<<{final_answer}>>>")


solve()