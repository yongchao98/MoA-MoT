import math

def calculate_integral():
    """
    This function calculates the definite integral by evaluating the antiderivative at the limits.
    """
    print("The definite integral is found by evaluating the antiderivative F(u) = 2^(17/16) * (sin(u))^(1/4) * tan(u) at the integration limits u = pi/4 and u = 0.")
    print("The calculation is F(pi/4) - F(0).\n")

    # --- Calculation for F(pi/4) ---
    print("--- Calculation for the upper limit F(pi/4) ---")
    
    # We want to compute: 2**(17/16) * (sin(pi/4))**(1/4) * tan(pi/4)
    c = "2^(17/16)"
    sin_pi4_str = "sin(pi/4)"
    tan_pi4_str = "tan(pi/4)"
    print(f"F(pi/4) = {c} * ({sin_pi4_str})^(1/4) * {tan_pi4_str}")

    sin_pi4_val = "2^(-1/2)"
    tan_pi4_val = "1"
    print(f"Since {sin_pi4_str} = 1/sqrt(2) = {sin_pi4_val} and {tan_pi4_str} = {tan_pi4_val}, we have:")
    
    print(f"F(pi/4) = 2^(17/16) * ({sin_pi4_val})^(1/4) * {tan_pi4_val}")

    exp_term_2 = "-1/8"
    print(f"         = 2^(17/16) * 2^({exp_term_2})")
    
    # Exponent arithmetic
    # 17/16 + (-1/8) = 17/16 - 2/16 = 15/16
    exp1_num, exp1_den = 17, 16
    exp2_num, exp2_den = -1, 8
    final_exp_num = 15
    final_exp_den = 16
    
    print(f"We add the exponents: {exp1_num}/{exp1_den} + ({exp2_num}/{exp2_den}) = 17/16 - 2/16 = {final_exp_num}/{final_exp_den}.")
    f_upper_str = f"2^({final_exp_num}/{final_exp_den})"
    print(f"So, F(pi/4) = {f_upper_str}\n")
    
    # --- Calculation for F(0) ---
    print("--- Calculation for the lower limit F(0) ---")
    f_lower_val = 0
    print(f"F(0) = 2^(17/16) * (sin(0))^(1/4) * tan(0) = {f_lower_val}\n")
    
    # --- Final Result ---
    print("--- Final Result ---")
    print(f"The value of the integral = F(pi/4) - F(0) = {f_upper_str} - {f_lower_val} = {f_upper_str}")
    
    # Numerical calculation
    final_answer = math.pow(2, 15/16)
    print(f"The numerical value is approximately: {final_answer:.7f}")

if __name__ == '__main__':
    calculate_integral()