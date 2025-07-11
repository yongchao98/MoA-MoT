import math

def calculate_shutdown_time():
    """
    Calculates the answers for all four parts of the problem.
    """
    # Given values
    P = 100  # Initial power in Watts
    x = 4    # Control subsystem power
    y = 5    # Camera subsystem power
    z = 1    # Sensor subsystem power
    T_half = 400 # Half-life in days

    # Part a: When to turn off the camera subsystem
    # Total power required with all systems on
    power_req_a = x + y + z
    # Calculate time using the derived formula: t = T_half * ln(P/P_req) / ln(2)
    t_a = T_half * math.log(P / power_req_a) / math.log(2)
    # Round to the nearest total day
    answer_a = round(t_a)

    # Part b: When to turn off the sensor subsystem
    # Total power required with camera off
    power_req_b = x + z
    # Calculate time
    t_b = T_half * math.log(P / power_req_b) / math.log(2)
    # Round to the nearest total day
    answer_b = round(t_b)

    # Part c: Memory usage in D
    # As determined by analyzing the hypothetical efficient C program:
    # 9 'frac' variables (6D each) + 2 'int' variables (5D each)
    # (9 * 6) + (2 * 5) = 54 + 10 = 64
    answer_c = 64

    # Part d: Number of exp function calls
    # As determined by analyzing the number of calls to ln(2), ln(10), ln(20)
    # calls = ceil(log2(2e9)) + ceil(log2(10e9)) + ceil(log2(20e9))
    # calls = 31 + 34 + 35 = 100
    answer_d = 100
    
    # Print the final answer in the format a:b:c:d
    print(f"a={answer_a}:b={answer_b}:c={answer_c}:d={answer_d}")

calculate_shutdown_time()
<<<1329:1729:64:100>>>