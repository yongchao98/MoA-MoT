import math

def solve_wuxing_probe_problem():
    """
    Solves the Pioneer probe power problem based on the Wuxing architecture.
    """
    # Input values for the specific case
    P = 100  # Initial power in Watts
    x = 4    # Control subsystem power in Watts
    y = 5    # Camera subsystem power in Watts
    z = 1    # Sensor subsystem power in Watts

    # --- Part a: Camera subsystem shutdown time ---
    # The camera is turned off when power drops below the total needed for all three systems.
    # We solve for t in: (x + y + z) = P * (1/2)^(t/400)
    # t = 400 * ln(P / (x + y + z)) / ln(2)
    p_target_camera = x + y + z
    t_camera_float = 400 * math.log(P / p_target_camera) / math.log(2)
    # Round to the nearest whole day
    a = int(t_camera_float + 0.5)

    # --- Part b: Sensor subsystem shutdown time ---
    # The sensor is turned off when power drops below what's needed for control and sensor.
    # We solve for t in: (x + z) = P * (1/2)^(t/400)
    # t = 400 * ln(P / (x + z)) / ln(2)
    p_target_sensor = x + z
    t_sensor_float = 400 * math.log(P / p_target_sensor) / math.log(2)
    # Round to the nearest whole day
    b = int(t_sensor_float + 0.5)

    # --- Part c: Memory usage in D (decimal digits) ---
    # We calculate the memory for variables in a memory-efficient C program.
    # Data type sizes: unsigned int = 5D, struct frac = 6D.
    # Variables needed:
    # 4 `unsigned int` for inputs P, x, y, z: 4 * 5D = 20D
    # 1 `unsigned int` for target power (reused): 1 * 5D = 5D
    # 1 `frac` to store the pre-calculated value of ln(2): 1 * 6D = 6D
    # 1 `frac` to store the result of the time calculation (reused): 1 * 6D = 6D
    c = (4 * 5) + (1 * 5) + (1 * 6) + (1 * 6)

    # --- Part d: Number of `exp` function calls ---
    # The ln(x) function calls exp(m) in a loop.
    # Number of iterations k for ln(x) is floor(log2(x * 1e9)) + 1.
    # The program needs to calculate ln(2), ln(P/p_target_camera)=ln(10), and ln(P/p_target_sensor)=ln(20).
    
    def count_exp_calls(val):
        """Calculates the number of exp calls needed to compute ln(val)."""
        if val <= 0:
            return 0
        # The number of iterations is determined by the binary search range and precision
        # (val / 2^k) < 1e-9  =>  val * 1e9 < 2^k  =>  log2(val * 1e9) < k
        return math.floor(math.log2(val * 1e9)) + 1

    calls_for_ln2 = count_exp_calls(2)
    calls_for_ln10 = count_exp_calls(P / p_target_camera) # ln(100/10) = ln(10)
    calls_for_ln20 = count_exp_calls(P / p_target_sensor) # ln(100/5) = ln(20)
    
    # Total calls in an efficient program that computes ln(2) once
    d = calls_for_ln2 + calls_for_ln10 + calls_for_ln20

    # The final equation string a:b:c:d
    final_answer = f"{a}:{b}:{c}:{d}"
    print(final_answer)

solve_wuxing_probe_problem()
<<<1329:1729:37:100>>>