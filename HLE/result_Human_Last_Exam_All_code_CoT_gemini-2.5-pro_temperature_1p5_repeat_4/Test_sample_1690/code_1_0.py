import math

def solve_pioneer_problem():
    """
    Solves the Pioneer probe problem by calculating shutdown times, memory usage,
    and function call counts, then prints the result in the specified format.
    """

    # --- Input Parameters ---
    P = 100.0       # Initial power in Watts
    x = 4.0         # Control subsystem power in Watts
    y = 5.0         # Camera subsystem power in Watts
    z = 1.0         # Sensor subsystem power in Watts
    T_half = 400.0  # Battery half-life in days

    # --- Part a: Camera Shutdown Time ---
    # Shutdown occurs when power drops to the level required by all systems.
    # P_target_a = x + y + z
    # The formula to find the time 't' is derived from P(t) = P * (0.5)^(t/T_half)
    # t = T_half * log2(P / P_target)
    P_target_a = x + y + z
    t_camera = T_half * math.log2(P / P_target_a)
    answer_a = round(t_camera)

    # --- Part b: Sensor Shutdown Time ---
    # Shutdown occurs when power drops to the level required by control and sensor systems.
    # P_target_b = x + z
    P_target_b = x + z
    t_sensor = T_half * math.log2(P / P_target_b)
    answer_b = round(t_sensor)

    # --- Part c: Memory Usage in D ---
    # Analysis of the most efficient Wuxing C program reveals 5 `frac` variables
    # are needed: P, a constant C = T_half/ln(2), a reusable P_target, a reusable
    # result 't', and T_half.
    # Size of `frac` = n(2D) + d(2D) + e(2D) = 6D.
    # Total memory = 5 variables * 6D/variable.
    answer_c = 30

    # --- Part d: Number of exp Function Calls ---
    # The efficient program calculates t = (T_half/ln(2)) * ln(P/P_target).
    # This requires computing ln(2), ln(P/P_target_a), and ln(P/P_target_b).
    # P/P_target_a = 100/10 = 10. P/P_target_b = 100/5 = 20.
    # We need to calculate ln(2), ln(10), and ln(20).
    # Number of exp calls for ln(v) is ceil(log2(v * 1e9)).
    calls_ln2 = math.ceil(math.log2(2 * 1e9))   # 31
    calls_ln10 = math.ceil(math.log2(10 * 1e9)) # 34
    calls_ln20 = math.ceil(math.log2(20 * 1e9)) # 35
    answer_d = calls_ln2 + calls_ln10 + calls_ln20 # 31 + 34 + 35 = 100

    # --- Final Result ---
    # Print the answers in the format a:b:c:d
    print(f"{answer_a}:{answer_b}:{answer_c}:{answer_d}")

solve_pioneer_problem()