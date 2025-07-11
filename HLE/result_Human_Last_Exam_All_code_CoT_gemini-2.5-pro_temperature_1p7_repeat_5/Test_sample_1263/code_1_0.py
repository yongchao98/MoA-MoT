import math

# Titan architecture constraint: numbers in fractions must be <= 15
MAX_VAL = 15

def check_overflow(num, den):
    """Checks if a fraction overflows the 4-bit constraint."""
    return num > MAX_VAL or den > MAX_VAL

def print_step(message, num, den, exp, success=True):
    """Helper function to print a calculation step."""
    status = "OK" if success else "OVERFLOW"
    print(f"{message:<45} | Expression: {num}/{den} * 10^{exp} ({status})")

def main():
    """
    Simulates the calculation of Pandora's escape velocity on the Titan architecture.
    """
    print("Objective: Calculate Escape Velocity v_e of Pandora.")
    print("Formula: v_e^2 = (8 * pi * G * R^2 * rho) / 3")
    print("We will simulate the calculation of the mantissa part: 8/3 * pi * R^2_mantissa * rho_mantissa\n")
    
    # --- Define Constants and Parameters as Titan Fractions ---
    # Radius R = 2000 km = 2e6 m. R^2 = 4e12. Mantissa is 4.
    r_squared = (4, 1, 12)
    # Shell density rho = 0.3 t/m^3 = 300 kg/m^3 = 3e2. Mantissa is 3.
    rho_shell = (3, 1, 2)
    # Gravitational constant G ~ 6.67e-11. Approx G = 13/2 e-11 = 6.5e-11
    G = (13, 2, -11)
    # Pi can be approximated as 3/1
    pi = (3, 1, 0)
    # The 8/3 term from the formula
    eight_thirds = (8, 3, 0)
    
    print("--- Titan Simulation ---")
    print("Using approximation pi = 3/1 to simplify the calculation.\n")

    # Let's try to compute M_mantissa = 4/3 * pi * R^3 * rho, which simplifies to 4 * R^3 * rho
    # M_mantissa requires R^3 mantissa (2e6 -> 2, so 2^3=8) and rho_shell mantissa (3).
    # Product is 4 * 8 * 3 = 96. Let's see if Titan can compute it step-by-step.
    
    # Start with the 4/3 factor in mass calculation M = 4/3 * pi * ...
    current_num, current_den, current_exp = 4, 3, 0
    print_step("MOV AX, 4/3", current_num, current_den, current_exp)

    # Multiply by pi = 3/1
    print_step(f"MUL AX, {pi[0]}/{pi[1]} (pi)", pi[0], pi[1], pi[2])
    # The intermediate product would be (4*3)/(3*1) = 12/3.
    # Titan's "algebraic simplification" capability allows cancelling common factors.
    # 12/3 simplifies to 4/1.
    current_num, current_den = 4, 1
    print_step("  -> RED AX (simplified)", current_num, current_den, current_exp, success=not check_overflow(current_num, current_den))
    
    # Now, multiply by the mantissa of R^3, which is 8
    r_cubed_mantissa = (8, 1, 18)
    print_step(f"MUL AX, {r_cubed_mantissa[0]}/{r_cubed_mantissa[1]} (R^3 mantissa)", r_cubed_mantissa[0], r_cubed_mantissa[1], 0)
    
    next_num = r_cubed_mantissa[0]
    result_num = current_num * next_num

    print(f"  -> Intermediate product: {current_num} * {next_num} = {result_num}")

    # Check for overflow
    if check_overflow(result_num, 1):
        print_step("  -> RED AX", result_num, 1, current_exp, success=False)
        print("\n--- ANALYSIS ---")
        print(f"The operation results in the numerator {result_num}, which is greater than {MAX_VAL}.")
        print("According to Titan rules, this must be simplified.")
        print("The provided simplification rule applies to fractions of the form B/C where B > C (e.g., 6/5 = 1 + 1/5).")
        print(f"The multiplication '{current_num}/1 * {next_num}/1' involves two integers. No simplification rule is defined for this case.")
        print("This overflow makes the calculation impossible on the Titan architecture.")
        print("\nFinal Conclusion: It is not possible to calculate the escape velocity.")
        
    else:
        # This part will not be reached, but is here for completeness
        current_num = result_num
        print_step("  -> RED AX", current_num, current_den, current_exp, success=True)

main()
print("\nFinal Answer:")
print("The final equation for the escape velocity squared is v_e^2 = 8/3 * pi * G * R^2 * rho_shell.")
print("To calculate this, we must multiply the contributing numbers: 8, 3, pi (e.g., 3), G (e.g., 13/2), R^2 (e.g., 4), rho (e.g., 3).")
print("For instance, a required intermediate step is 8 * 4 = 32.")
print("Since 32 is greater than 15, it overflows the 4-bit integer limit.")
print("The rules provided for Titan do not specify a method for simplifying the product of two integers that results in an overflow.")
print("Therefore, the calculation is impossible.")