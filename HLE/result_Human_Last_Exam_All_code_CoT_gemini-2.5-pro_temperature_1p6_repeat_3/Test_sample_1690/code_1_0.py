import math

def solve_pioneer_problem():
    """
    Solves the Pioneer power decay problem and analyzes the
    corresponding efficient Wuxing C program.
    """
    # Initial given values
    P = 100  # Initial power in Watts
    x = 4    # Control subsystem power in Watts
    y = 5    # Camera subsystem power in Watts
    z = 1    # Sensor subsystem power in Watts
    T_half = 400 # Power half-life in days

    # --- Part a: Time to turn off the camera subsystem ---
    # The camera turns off when the available power is no longer sufficient for all three subsystems.
    # The equation to solve is: P * (1/2)^(t / T_half) = x + y + z
    # Solving for t: t = T_half * ln(P / (x + y + z)) / ln(2)
    P_req_camera = x + y + z
    t_camera_float = T_half * math.log(P / P_req_camera) / math.log(2)
    a = round(t_camera_float)

    # --- Part b: Time to turn off the sensor subsystem ---
    # The sensor turns off when power is no longer sufficient for the control and sensor subsystems.
    # The equation to solve is: P * (1/2)^(t / T_half) = x + z
    # Solving for t: t = T_half * ln(P / (x + z)) / ln(2)
    P_req_sensor = x + z
    t_sensor_float = T_half * math.log(P / P_req_sensor) / math.log(2)
    b = round(t_sensor_float)

    # --- Part c: Memory usage analysis ---
    # For the most time-and-memory-efficient C program, we need to minimize variables
    # while also reducing expensive calculations. An optimal approach involves reusing calculated
    # logarithm values.
    # Variables in main():
    # 5 `frac` inputs (P, x, y, z, T_half)
    # 1 `frac` for a temporary result variable
    # 1 `frac` to store the result of ln(2)
    # 1 `frac` to store the result of ln(10)
    # Total `frac` variables = 8.
    # Memory for `frac` type = 3 * sizeof(char) = 3 * 2D = 6D.
    # Total memory usage = 8 variables * 6D/variable = 48D.
    c = 48

    # --- Part d: `exp` function call analysis ---
    # The most time-efficient program calculates ln(2) and ln(10), then derives ln(20)
    # from ln(2) + ln(10) to avoid a third expensive ln() call.
    # The number of `exp` calls inside ln(N) is ceil(log2(N * 10^9)).
    # 1. Calls for ln(2): ceil(log2(2 * 10^9)) = 31 calls.
    # 2. Calls for ln(10): ceil(log2(10 * 10^9)) = 34 calls.
    # Total `exp` calls = 31 + 34 = 65.
    d = 65

    # Print the final result in the format a:b:c:d
    print(f"a. Shutdown time for camera subsystem (rounded to days):")
    print(f"   Equation: t = {T_half} * ln({P} / ({x} + {y} + {z})) / ln(2)")
    print(f"   Result: {a}\n")

    print(f"b. Shutdown time for sensor subsystem (rounded to days):")
    print(f"   Equation: t = {T_half} * ln({P} / ({x} + {z})) / ln(2)")
    print(f"   Result: {b}\n")
    
    print(f"c. Memory usage for variables in the efficient Wuxing program:")
    print(f"   Analysis: 8 frac variables * 6D/frac")
    print(f"   Result: {c}\n")

    print(f"d. Number of 'exp' calls in the efficient Wuxing program:")
    print(f"   Analysis: Calls for ln(2) + Calls for ln(10) = 31 + 34")
    print(f"   Result: {d}\n")
    
    print("---")
    print(f"Final answer in a:b:c:d format: {a}:{b}:{c}:{d}")


solve_pioneer_problem()