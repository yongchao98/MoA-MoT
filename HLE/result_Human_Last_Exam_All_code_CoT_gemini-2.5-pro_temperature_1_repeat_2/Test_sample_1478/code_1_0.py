import math

# Titan architecture constants and helper functions

# The maximum value for a 6-bit unsigned integer is 2^6 - 1 = 63.
MAX_VAL = 63

def check_overflow(val):
    """Checks if a value exceeds the 6-bit limit."""
    if val > MAX_VAL:
        return True
    return False

def simplify_fraction(n, d):
    """Simplifies a fraction by dividing by the greatest common divisor."""
    if d == 0:
        raise ZeroDivisionError("Denominator cannot be zero.")
    common_divisor = math.gcd(n, d)
    return n // common_divisor, d // common_divisor

def print_step(instruction, register, value, comment=""):
    """Formats and prints a single step of the Titan program."""
    print(f"{instruction:<10} {register:<5}, {value:<10} ; {comment}")

def print_register_state(register, expression_str, comment=""):
    """Prints the current state of a register."""
    print(f"; {register} holds the expression: {expression_str}")
    if comment:
        print(f"; {comment}")
    print("-" * 50)

def main():
    """
    Simulates the calculation of the gravitational force on the Titan architecture.
    """
    print("Starting Titan program to calculate gravitational force.")
    print("All intermediate numerators and denominators must be <= 63.\n")

    # F = G * M * m / d^2
    # M = rho * V = rho * (4/3) * pi * R^3
    # F = G * m * rho * (4/3) * pi * R^3 / d^2
    # We will calculate the product of the mantissas.
    # The exponent calculation is separate:
    # G_exp = -11, m_exp = 0, rho_exp = 3, R^3_exp = 18, d^-2_exp = -6
    # Total_exp = -11 + 3 + 18 - 6 = 4.

    # Mantissas of the constants:
    # G: 20/3
    # m: 50/1
    # rho: 6/5
    # 4/3: 4/3
    # pi: 22/7
    # R^3: 8/1
    # d^2: 1/1 (exponent handles this)

    print("Plan: Multiply mantissas sequentially and check for overflow at each step.")
    print("Let's compute the product of (4/3) * R^3_m * rho_m * G_m first.")
    print("-" * 50)

    # Step 1: Start with 4/3
    n, d = 4, 3
    print_step("MOV", "AX", "4/3", "Load volumetric constant")
    print_register_state("AX", f"{n}/{d}")
    
    # Step 2: Multiply by R^3 mantissa (8)
    mul_n, mul_d = 8, 1
    print_step("MUL", "AX", f"{mul_n}/{mul_d}", "Multiply by R^3 mantissa (8)")
    n, d = n * mul_n, d * mul_d
    n, d = simplify_fraction(n, d)
    
    if check_overflow(n) or check_overflow(d):
        print(f"; Operation failed: Result {n}/{d} contains a value > {MAX_VAL}")
        print("<<<N0>>>")
        return
        
    print_register_state("AX", f"{n}/{d}", "Current value of (4/3)*8 is 32/3.")
    
    # Step 3: Multiply by rho mantissa (6/5)
    mul_n, mul_d = 6, 5
    print_step("MUL", "AX", f"{mul_n}/{mul_d}", "Multiply by density mantissa (6/5)")
    n, d = n * mul_n, d * mul_d
    n, d = simplify_fraction(n, d)

    if check_overflow(n) or check_overflow(d):
        print(f"; Operation failed: Result {n}/{d} contains a value > {MAX_VAL}")
        print("<<<N0>>>")
        return

    print_register_state("AX", f"{n}/{d}", "Current value is 64/5.")

    print("; Here, the calculation fails.")
    print(f"; The resulting numerator is {n}, which is not representable as a 6-bit integer (0-63).")
    print("; According to Titan rules, this operation is invalid.")
    print("; We cannot use decomposition strategies like (32/3) * (1 + 1/5) because the intermediate product")
    print("; of the integer parts (32 * 1) is fine, but the fractional part (32/3)*(1/5) = 32/15 is added back,")
    print("; and the next multiplications will also overflow.")
    
    print("\n; Let's demonstrate the simplest case of this failure.")
    print("; The calculation requires multiplying 32 by other factors, like 2.")
    print("-" * 50)
    
    n_fail, d_fail = 32, 1
    print_step("MOV", "CX", f"{n_fail}/{d_fail}", "Isolate the problematic component")
    print_register_state("CX", f"{n_fail}/{d_fail}")

    mul_n_fail, mul_d_fail = 2, 1
    print_step("MUL", "CX", f"{mul_n_fail}/{mul_d_fail}", "Multiply by smallest possible factor > 1")
    
    res_n = n_fail * mul_n_fail
    
    print(f"; Attempting to compute {n_fail} * {mul_n_fail} = {res_n}")
    print(f"; The result {res_n} exceeds the max 6-bit value of {MAX_VAL}.")
    print("; Since 2 is prime, this operation cannot be decomposed further.")
    print("; The Titan architecture cannot perform this calculation with the given initial values.")
    
    print("\nConclusion: The calculation is not possible.")
    print("Final Answer: N0")
    print("<<<N0>>>")


if __name__ == "__main__":
    main()
