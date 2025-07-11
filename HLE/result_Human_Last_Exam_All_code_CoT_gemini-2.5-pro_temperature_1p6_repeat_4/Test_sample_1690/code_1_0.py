import math

def solve_pioneer_problem():
    """
    Solves the Pioneer probe power problem based on the Wuxing architecture specifications.
    """
    # Given inputs
    P = 100  # Initial power (W)
    x = 4    # Control subsystem power (W)
    y = 5    # Camera subsystem power (W)
    z = 1    # Sensor subsystem power (W)
    T_half = 400 # Power half-life (days)

    print("--- Solving for given values: P=100, x=4, y=5, z=1 ---")

    # a. When to turn off the camera subsystem?
    print("\n[Part a: Camera Shutdown Time]")
    p_final_cam = x + y + z
    ratio_cam = P / p_final_cam
    
    # The equation to solve is: t = T_half * ln(P / (x+y+z)) / ln(2)
    print(f"Equation: t_cam = {T_half} * ln({P} / ({x}+{y}+{z})) / ln(2)")
    print(f"Equation: t_cam = {T_half} * ln({ratio_cam}) / ln(2)")

    t_cam = T_half * math.log(ratio_cam) / math.log(2)
    a_answer = round(t_cam)
    print(f"Result: t_cam = {t_cam:.2f} days, which rounds to {a_answer} days.")

    # b. When to turn off the sensor subsystem?
    print("\n[Part b: Sensor Shutdown Time]")
    p_final_sensor = x + z
    ratio_sensor = P / p_final_sensor
    
    # The equation to solve is: t = T_half * ln(P / (x+z)) / ln(2)
    print(f"Equation: t_sensor = {T_half} * ln({P} / ({x}+{z})) / ln(2)")
    print(f"Equation: t_sensor = {T_half} * ln({ratio_sensor}) / ln(2)")
    
    t_sensor = T_half * math.log(ratio_sensor) / math.log(2)
    b_answer = round(t_sensor)
    print(f"Result: t_sensor = {t_sensor:.2f} days, which rounds to {b_answer} days.")
    
    # c. What is the memory usage in D for variables of this program?
    print("\n[Part c: Memory Usage Analysis]")
    # Wuxing data type sizes: int=5D, frac=6D (3*2D)
    # Minimal set of variables for a memory-efficient C program:
    # - Inputs: int P, x, y, z; (4 ints)
    # - Results: int t_cam, t_sensor; (2 ints)
    # - Intermediate values for calculation:
    #   - int p_final; (1 int, can be reused)
    #   - frac ln2, ratio, t_frac; (3 fracs, 'ratio' and 't_frac' are reused)
    # Total ints: 4 (inputs) + 2 (results) + 1 (intermediate) = 7 ints
    # Total fracs: 3 (intermediate)
    int_mem = 7 * 5
    frac_mem = 3 * 6
    c_answer = int_mem + frac_mem
    print("Variable memory allocation:")
    print(f"- 7 int variables (inputs, results, intermediate): 7 * 5D = {int_mem}D")
    print(f"- 3 frac variables (ln2, ratio, t_frac): 3 * 6D = {frac_mem}D")
    print(f"Total memory usage = {int_mem}D + {frac_mem}D = {c_answer}D")

    # d. What is the number of times this program calls function exp?
    print("\n[Part d: exp() Call Analysis]")
    # exp() is called once per iteration of the ln(v) function's binary search.
    # The ln(v) loop runs k times, where k = ceil(log2(v) + log2(10^9)).
    # log2(10^9) = 9 * log2(10) approx 9 * 3.3219 = 29.897
    # So k = ceil(log2(v) + 29.897)
    # The program calculates ln(2), ln(ratio_cam) which is ln(10), and ln(ratio_sensor) which is ln(20).
    log2_10_9 = 9 * math.log2(10)

    # For ln(2)
    k_ln2 = math.ceil(math.log2(2) + log2_10_9)
    # For ln(10)
    k_ln10 = math.ceil(math.log2(10) + log2_10_9)
    # For ln(20)
    k_ln20 = math.ceil(math.log2(20) + log2_10_9)
    
    d_answer = k_ln2 + k_ln10 + k_ln20
    print(f"exp() calls for ln(2): {k_ln2}")
    print(f"exp() calls for ln(10): {k_ln10}")
    print(f"exp() calls for ln(20): {k_ln20}")
    print(f"Total exp() calls = {k_ln2} + {k_ln10} + {k_ln20} = {d_answer}")

    # Final combined answer
    print("\n--- Final Answer ---")
    final_answer_string = f"{a_answer}:{b_answer}:{c_answer}:{d_answer}"
    print(final_answer_string)
    
    # Return for format requirement
    return f"<<<{final_answer_string}>>>"

# Execute the function and print the final formatted answer
final_answer = solve_pioneer_problem()
# The final formatted answer is part of the output above, but we also capture it
# to present at the very end as requested.
# print(final_answer) # This would print it again. The structure is self-contained.
