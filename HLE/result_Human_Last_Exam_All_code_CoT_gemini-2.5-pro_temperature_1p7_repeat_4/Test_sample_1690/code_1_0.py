import math

def solve():
    """
    Solves the Pioneer probe problem based on the Wuxing architecture specs.
    """
    # Input values
    P = 100
    x = 4
    y = 5
    z = 1
    T_half = 400

    # --- Part a & b: Calculate Shutdown Times ---
    
    # a. Camera shutdown time
    power_threshold_a = x + y + z
    ratio_a = P / power_threshold_a
    # Formula: t_a = T_half * ln(ratio_a) / ln(2)
    t_a_float = T_half * math.log(ratio_a) / math.log(2)
    t_a = round(t_a_float)

    print("--- Calculations ---")
    print("\na. Time to turn off the camera subsystem:")
    print(f"Power threshold = x + y + z = {x} + {y} + {z} = {power_threshold_a} W")
    print(f"t_a = {T_half} * ln({P} / {power_threshold_a}) / ln(2) = {T_half} * ln({ratio_a}) / ln(2)")
    print(f"t_a ≈ {t_a_float:.2f} days, which is rounded to {t_a} days.")

    # b. Sensor shutdown time
    power_threshold_b = x + z
    ratio_b = P / power_threshold_b
    # Formula: t_b = T_half * ln(ratio_b) / ln(2)
    t_b_float = T_half * math.log(ratio_b) / math.log(2)
    t_b = round(t_b_float)
    
    print("\nb. Time to turn off the sensor subsystem:")
    print(f"Power threshold = x + z = {x} + {z} = {power_threshold_b} W")
    print(f"t_b = {T_half} * ln({P} / {power_threshold_b}) / ln(2) = {T_half} * ln({ratio_b}) / ln(2)")
    print(f"t_b ≈ {t_b_float:.2f} days, which is rounded to {t_b} days.")

    # --- Part c: Calculate Memory Usage ---
    
    # Minimal variables for an efficient program:
    # 4 for inputs (P, x, y, z)
    # 1 for ln(2) result (ln2_val)
    # 1 for ln(10) result (ln10_val)
    # 1 for the final calculated time (result_t)
    # Total = 7 variables
    num_frac_vars = 7
    # Size of one frac = size(n) + size(d) + size(e) = 2D + 2D + 2D = 6D
    mem_per_frac = 6
    mem_usage = num_frac_vars * mem_per_frac

    print("\nc. Memory usage in D for variables of this program:")
    print(f"The most efficient program stores inputs (P,x,y,z), intermediate results (ln(2), ln(10)) and a final result variable.")
    print(f"Total variables needed = {num_frac_vars} (all of type frac)")
    print(f"Size of one 'frac' variable = 2D (n) + 2D (d) + 2D (e) = {mem_per_frac}D")
    print(f"Total memory usage = {num_frac_vars} variables * {mem_per_frac} D/variable = {mem_usage}D")

    # --- Part d: Calculate exp() calls ---

    # An efficient program calls ln(2) and ln(10).
    # Calls for ln(x) = ceil(log2(x * 1e9))
    exp_calls_for_ln2 = math.ceil(math.log2(2 * 1e9))
    exp_calls_for_ln10 = math.ceil(math.log2(10 * 1e9))
    total_exp_calls = exp_calls_for_ln2 + exp_calls_for_ln10
    
    print("\nd. Number of times this program calls function exp:")
    print("An efficient program calls ln(2) and ln(10) once.")
    print(f"Number of exp() calls for ln(2) ≈ ceil(log2(2 * 10^9)) = {exp_calls_for_ln2}")
    print(f"Number of exp() calls for ln(10) ≈ ceil(log2(10 * 10^9)) = {exp_calls_for_ln10}")
    print(f"Total exp() calls = {exp_calls_for_ln2} + {exp_calls_for_ln10} = {total_exp_calls}")
    
    # --- Final Answer ---
    final_answer = f"{t_a}:{t_b}:{mem_usage}:{total_exp_calls}"
    print("\n--- Final Answer ---")
    print(f"The final answer in the format a:b:c:d is:\n\n<<<{final_answer}>>>")

solve()