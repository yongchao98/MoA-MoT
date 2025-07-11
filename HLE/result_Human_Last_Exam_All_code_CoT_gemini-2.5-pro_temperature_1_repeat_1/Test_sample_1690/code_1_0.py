import math

def solve_pioneer_problem():
    """
    Calculates the answers to the Pioneer probe problem based on the Wuxing architecture.
    """
    # Given values
    P = 100  # Initial power in Watts
    x = 4    # Control subsystem power in Watts
    y = 5    # Camera subsystem power in Watts
    z = 1    # Sensor subsystem power in Watts
    T_half = 400 # Half-life in days

    # Part a: Calculate when to turn off the camera subsystem
    # This happens when power drops below the requirement for all three systems.
    power_threshold_a = x + y + z
    # The formula is t = T_half * ln(P_initial / P_threshold) / ln(2)
    # t_a = 400 * ln(100 / 10) / ln(2) = 400 * ln(10) / ln(2)
    t_a_float = T_half * math.log(P / power_threshold_a) / math.log(2)
    a = round(t_a_float)

    # Part b: Calculate when to turn off the sensor subsystem
    # This happens when power drops below the requirement for control and sensor systems.
    power_threshold_b = x + z
    # t_b = 400 * ln(100 / 5) / ln(2) = 400 * ln(20) / ln(2)
    t_b_float = T_half * math.log(P / power_threshold_b) / math.log(2)
    b = round(t_b_float)

    # Part c: Calculate the memory usage in D for variables of the program
    # A time-and-memory efficient C program would compute ln(2) once and reuse it.
    # The program would calculate t1 and t2.
    # Minimal variables needed in main(): frac ln2, frac t1, frac t2. Total = 3 variables.
    # Size of one 'frac' variable = sizeof(n) + sizeof(d) + sizeof(e)
    # sizeof(char) on Wuxing = 2D
    # Size of frac = 2D + 2D + 2D = 6D.
    # Total memory = 3 variables * 6D/variable.
    c = 3 * 6

    # Part d: Calculate the number of times the program calls function exp
    # The most efficient program calls ln(2), ln(10), and ln(20).
    # The ln(x) function calls exp() in a loop. Number of calls = ceil(log2(x / 1e-9)).
    epsilon = 1e-9
    # Calls for ln(2)
    exp_calls_ln2 = math.ceil(math.log2(2 / epsilon))
    # Calls for ln(10) (used for t_a calculation)
    exp_calls_ln10 = math.ceil(math.log2(power_threshold_a / epsilon))
    # Calls for ln(20) (used for t_b calculation)
    exp_calls_ln20 = math.ceil(math.log2(P / power_threshold_b / epsilon))
    
    # Total exp() calls is the sum of calls from each ln() invocation.
    d = exp_calls_ln2 + exp_calls_ln10 + exp_calls_ln20

    # Print the final answer in the format a:b:c:d
    print(f"{a}:{b}:{c}:{d}")

solve_pioneer_problem()