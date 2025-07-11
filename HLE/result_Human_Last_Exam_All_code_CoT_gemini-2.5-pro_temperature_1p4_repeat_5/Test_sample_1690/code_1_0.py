import math

def solve_pioneer_problem():
    """
    Calculates the shutdown times and program metrics for the Pioneer probe
    on the Wuxing architecture.
    """
    # --- Givens ---
    P = 100.0
    x = 4.0
    y = 5.0
    z = 1.0
    half_life_days = 400.0

    print("--- Solving for Pioneer Probe ---")
    print(f"Initial Power (P): {P}W")
    print(f"Subsystem Power -> Control (x): {x}W, Camera (y): {y}W, Sensor (z): {z}W")
    print(f"Power Half-life: {half_life_days} days\n")

    # --- Part a: Camera Shutdown Time ---
    # The shutdown happens when power drops below the total needed for all systems.
    p_camera_threshold = x + y + z
    # The formula is t = half_life * log2(P_initial / P_threshold)
    t_camera_exact = half_life_days * math.log2(P / p_camera_threshold)
    a = round(t_camera_exact)

    print("a. Camera Shutdown Time")
    print(f"   Shutdown Threshold: x + y + z = {x} + {y} + {z} = {p_camera_threshold} W")
    print(f"   Equation: t = {half_life_days} * log2({P} / {p_camera_threshold})")
    print(f"   Result: {t_camera_exact:.2f} days, rounded to {a} days.\n")

    # --- Part b: Sensor Shutdown Time ---
    # The shutdown happens when power drops below the total for control and sensor.
    p_sensor_threshold = x + z
    t_sensor_exact = half_life_days * math.log2(P / p_sensor_threshold)
    b = round(t_sensor_exact)

    print("b. Sensor Shutdown Time")
    print(f"   Shutdown Threshold: x + z = {x} + {z} = {p_sensor_threshold} W")
    print(f"   Equation: t = {half_life_days} * log2({P} / {p_sensor_threshold})")
    print(f"   Result: {t_sensor_exact:.2f} days, rounded to {b} days.\n")

    # --- Part c: Memory Usage in D ---
    # In an efficient Wuxing C program, we minimize variable declarations.
    # Sizes: frac = 6D, int = 5D
    # Vars: frac P,x,y,z (4*6D); int t_cam,t_sen (2*5D); frac ln2,lnP (2*6D);
    #       frac p_thresh, ln_val, t_frac (3*6D reused vars)
    mem_inputs = 4 * 6
    mem_outputs = 2 * 5
    mem_cached_calcs = 2 * 6
    mem_reused_vars = 3 * 6
    c = mem_inputs + mem_outputs + mem_cached_calcs + mem_reused_vars

    print("c. Memory Usage in D (Decimal Digits)")
    print("   Variable Memory Breakdown:")
    print(f"   - Inputs (4 frac): {mem_inputs}D | Outputs (2 int): {mem_outputs}D")
    print(f"   - Cached Calcs (2 frac): {mem_cached_calcs}D | Reused Vars (3 frac): {mem_reused_vars}D")
    print(f"   Equation: {mem_inputs} + {mem_outputs} + {mem_cached_calcs} + {mem_reused_vars}")
    print(f"   Result: {c} D\n")


    # --- Part d: Number of exp() Calls ---
    # The ln(X) function calls exp() k = ceil(log2(X * 10^9)) times.
    # Efficient code calls ln() for 2, 100, 10, and 5.
    def count_exp_calls(val):
        return math.ceil(math.log2(val * 1e9))

    calls_for_ln2 = count_exp_calls(2)
    calls_for_ln100 = count_exp_calls(P)
    calls_for_ln10 = count_exp_calls(p_camera_threshold)
    calls_for_ln5 = count_exp_calls(p_sensor_threshold)
    d = calls_for_ln2 + calls_for_ln100 + calls_for_ln10 + calls_for_ln5

    print("d. Total 'exp' Function Calls")
    print("   Calls from ln(2) = {}".format(calls_for_ln2))
    print("   Calls from ln(P=100) = {}".format(calls_for_ln100))
    print("   Calls from ln(camera_threshold=10) = {}".format(calls_for_ln10))
    print("   Calls from ln(sensor_threshold=5) = {}".format(calls_for_ln5))
    print(f"   Equation: {calls_for_ln2} + {calls_for_ln100} + {calls_for_ln10} + {calls_for_ln5}")
    print(f"   Result: {d} calls\n")

    # --- Final Combined Answer ---
    print("--- Final Result (a:b:c:d) ---")
    print(f"{a}:{b}:{c}:{d}")


if __name__ == '__main__':
    solve_pioneer_problem()