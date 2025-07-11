import math

def present_calculation():
    """
    Prints a step-by-step calculation that can be performed using integers and fractions.
    """
    # The constants for the problem
    d_int = 240
    g_int_num = 10
    tan_int_num = 27
    tan_int_den = 64

    # Determine the largest integer for the final answer
    largest_integer = max(d_int, g_int_num, tan_int_num, tan_int_den)
    answer_code = f"Y{largest_integer}"

    # --- Present the calculation for the son ---
    print("Yes, your son can calculate the time. Here is a simple way to frame the calculation for him:\n")
    print("To find the time `t`, we first need to calculate a value, let's call it `t-squared`, using this equation:")
    print(f"t-squared = (2 * {d_int} / {g_int_num}) * ({tan_int_num} / {tan_int_den})")
    print("\nLet's solve this step by step:\n")

    print("1. First, let's calculate the part in the first parenthesis:")
    calc1_res = 2 * d_int // g_int_num
    print(f"   2 * {d_int} / {g_int_num} = {2 * d_int} / {g_int_num} = {calc1_res}")

    print("\n2. Now our equation is much simpler:")
    print(f"   t-squared = {calc1_res} * ({tan_int_num} / {tan_int_den})")

    print("\n3. To work with smaller numbers, let's simplify by finding common factors.")
    print(f"   We can simplify '{calc1_res}' and '{tan_int_den}' by dividing both by 16:")
    print(f"   {calc1_res} / 16 = {calc1_res // 16}")
    print(f"   {tan_int_den} / 16 = {tan_int_den // 16}")
    calc2_num = calc1_res // 16
    calc2_den = tan_int_den // 16

    print("\n4. Now the multiplication is very easy:")
    t_sq_num = calc2_num * tan_int_num
    t_sq_den = calc2_den
    print(f"   t-squared = ({calc2_num} * {tan_int_num}) / {calc2_den} = {t_sq_num} / {t_sq_den}")

    print(f"\n5. We've found that t-squared = {t_sq_num} / {t_sq_den}. To find the time `t`, we just take the square root of the top and bottom:")
    t_num = int(math.sqrt(t_sq_num))
    t_den = int(math.sqrt(t_sq_den))
    print(f"   t = sqrt({t_sq_num}) / sqrt({t_sq_den}) = {t_num} / {t_den}")

    print(f"\nSo, the final answer is {t_num}/{t_den} seconds, which is 4.5 seconds.")
    
    print(f"\n<<<Y{largest_integer}>>>")

present_calculation()