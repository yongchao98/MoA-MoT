import math

def solve_pioneer_problem():
    """
    Calculates the shutdown times for the Pioneer probe's subsystems
    and provides analysis on the resource usage of a hypothetical Wuxing C program.
    """
    # Input parameters
    P = 100.0  # Initial power in Watts
    x = 4.0    # Control subsystem power in Watts
    y = 5.0    # Camera subsystem power in Watts
    z = 1.0    # Sensor subsystem power in Watts
    
    # Constants
    T_half = 400.0 # Half-life in days

    # a. Calculate when to turn off the camera subsystem
    # This happens when power drops to the level required for all three systems.
    power_target_cam = x + y + z
    # Using the formula: t = T_half * log2(P_initial / P_target)
    t_cam = T_half * math.log2(P / power_target_cam)
    # Round to the nearest whole day
    t_cam_rounded = int(t_cam + 0.5)

    # b. Calculate when to turn off the sensor subsystem
    # This happens when power drops to the level required for control and sensor.
    power_target_sensor = x + z
    t_sensor = T_half * math.log2(P / power_target_sensor)
    # Round to the nearest whole day
    t_sensor_rounded = int(t_sensor + 0.5)

    # c. Memory usage in D for variables of this program
    # Analysis:
    # An efficient program would have a main function and call ln(), which calls exp().
    # Peak variable count on stack:
    # - main: 10 frac vars (inputs, outputs, constants, temps)
    # - ln(): 4 frac vars (l, r, m, y)
    # - exp(): 2 frac vars (y, p), 1 char var (n)
    # Total peak: 16 frac, 1 char
    # Sizes: frac = 6D, char = 2D
    # Memory = 16 * 6D + 1 * 2D = 98D
    memory_usage_d = 98

    # d. Number of times this program calls function exp
    # Analysis:
    # An efficient program computes t = 400 * ln(P/P_target) / ln(2).
    # This requires 3 ln calls: ln(2), ln(P/(x+y+z))=ln(10), ln(P/(x+z))=ln(20).
    # Number of exp calls for ln(x) is k ≈ log2(x * 1e9).
    # k_ln2   = ceil(log2(2e9))  = 31
    # k_ln10  = ceil(log2(10e9)) = 34
    # k_ln20  = ceil(log2(20e9)) = 35
    # Total exp calls = 31 + 34 + 35 = 100
    exp_call_count = 100

    # Print the final answer in the required format a:b:c:d
    print(f"{t_cam_rounded}:{t_sensor_rounded}:{memory_usage_d}:{exp_call_count}")

solve_pioneer_problem()
<<<1329:1729:98:100>>>