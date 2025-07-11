import math

def solve_pioneer_problem():
    """
    Calculates the shutdown times for the Pioneer probe's subsystems
    and analyzes the theoretical Wuxing C program.
    """
    # Given parameters
    P = 100.0  # Initial power in Watts
    x = 4.0    # Control subsystem power in Watts
    y = 5.0    # Camera subsystem power in Watts
    z = 1.0    # Sensor subsystem power in Watts
    half_life = 400.0  # Power half-life in days

    # a. Camera subsystem shutdown time
    power_threshold_a = x + y + z
    # Equation: power_threshold_a = P * (0.5)^(t_a / half_life)
    # Solved for t_a: t_a = half_life * ln(power_threshold_a / P) / ln(0.5)
    t_a = half_life * math.log(power_threshold_a / P) / math.log(0.5)
    rounded_t_a = round(t_a)

    print(f"a: t = {int(half_life)} * ln({int(power_threshold_a)}/{int(P)}) / ln(0.5) = {rounded_t_a}")

    # b. Sensor subsystem shutdown time
    power_threshold_b = x + z
    # Equation: power_threshold_b = P * (0.5)^(t_b / half_life)
    # Solved for t_b: t_b = half_life * ln(power_threshold_b / P) / ln(0.5)
    t_b = half_life * math.log(power_threshold_b / P) / math.log(0.5)
    rounded_t_b = round(t_b)
    
    print(f"b: t = {int(half_life)} * ln({int(power_threshold_b)}/{int(P)}) / ln(0.5) = {rounded_t_b}")

    # c. Memory usage in D for variables of the Wuxing C program
    # A 'frac' variable is 3 * 'char' = 3 * 2D = 6D.
    # An efficient program uses:
    # 4 frac for inputs (P,x,y,z) = 24D
    # 4 frac for calculations (reused) = 24D
    # 1 frac for ln(0.5) constant = 6D
    # 1 temporary frac for '1/2e0' argument = 6D
    # Total = 24 + 24 + 6 + 6 = 60D
    memory_usage_c = 60
    
    print(f"c: Memory usage = {memory_usage_c}D")

    # d. Number of times the program calls function exp
    # The program calls ln() 3 times: ln(0.5), ln(0.1), and ln(0.05).
    # exp() calls for ln(x) = ceil(log2(x / 1e-9))
    exp_calls_ln0_5 = math.ceil(math.log2(0.5 / 1e-9))  # ~28.9 -> 29
    exp_calls_ln0_1 = math.ceil(math.log2(0.1 / 1e-9))  # ~26.57 -> 27
    exp_calls_ln0_05 = math.ceil(math.log2(0.05 / 1e-9)) # ~25.57 -> 26
    total_exp_calls_d = exp_calls_ln0_5 + exp_calls_ln0_1 + exp_calls_ln0_05
    
    print(f"d: Number of exp calls = {total_exp_calls_d}")

    # Print the final combined answer
    final_answer = f"{rounded_t_a}:{rounded_t_b}:{memory_usage_c}:{total_exp_calls_d}"
    print(f"\nFinal Answer: {final_answer}")

if __name__ == '__main__':
    solve_pioneer_problem()