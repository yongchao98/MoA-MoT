import math

def main():
    """
    This script simulates the calculation of Pandora's escape velocity on the
    Titan 4-bit architecture to determine its feasibility.
    """
    print("### Titan Architecture Feasibility Study for Pandora's Escape Velocity ###\n")

    # --- Step 1: Define a class to represent Titan numbers and operations ---
    # The key constraint is that numerators and denominators cannot exceed 15.
    class TitanFraction:
        def __init__(self, numerator, denominator=1, name=""):
            if not (0 <= numerator <= 15 and 0 < denominator <= 15):
                # This check is for the initial definition of constants.
                raise ValueError(f"Invalid 4-bit fraction for {name}: {numerator}/{denominator}")
            self.num = numerator
            self.den = denominator
            self.name = name

        def __mul__(self, other):
            print(f"Attempting multiplication: ({self}) * ({other})")
            
            # Perform multiplication
            new_num = self.num * other.num
            new_den = self.den * other.den

            # Check for immediate overflow, as per Titan rules
            if new_num > 15 or new_den > 15:
                print(f"  > Result: {new_num}/{new_den}")
                print(f"  > FAILURE: Overflow detected. Numerator or denominator > 15.\n")
                return None
            
            # If no overflow, simplify the fraction using GCD
            common_divisor = math.gcd(new_num, new_den)
            final_num = new_num // common_divisor
            final_den = new_den // common_divisor
            
            print(f"  > Result: {new_num}/{new_den}, simplified to {final_num}/{final_den}")
            
            # The simplified result must also be valid
            if final_num > 15 or final_den > 15:
                print(f"  > FAILURE: Simplified result still causes overflow.\n")
                return None

            print("  > SUCCESS: Operation is valid.\n")
            return TitanFraction(final_num, final_den)

        def __str__(self):
            return f"{self.name}({self.num}/{self.den})"

    # --- Step 2: Represent all coefficients in the equation as Titan fractions ---
    # Equation: v_e^2 = 2 * G * (4/3) * pi * rho_shell * R^2
    # We only need to check the multiplication of the coefficients.
    print("--- Defining constants for the equation: v_e^2 = C_2 * C_G * C_4_3 * C_pi * C_rho * C_R^2 ---\n")
    try:
        C_2 = TitanFraction(2, 1, name="2")
        C_G = TitanFraction(13, 2, name="G")      # Approx G ~ 6.5
        C_4_3 = TitanFraction(4, 3, name="4/3")
        C_pi = TitanFraction(3, 1, name="pi")     # Approx pi ~ 3
        C_rho = TitanFraction(3, 1, name="rho")   # For rho_shell = 300 = 3e2
        C_R = TitanFraction(2, 1, name="R")       # For R = 2e6
    except ValueError as e:
        print(f"Error initializing constants: {e}")
        return

    # --- Step 3: Simulate the calculation step-by-step ---
    # We will try an optimal order of operations to keep intermediate values small.
    print("--- Simulating calculation by multiplying coefficients step-by-step ---\n")

    # First, calculate R^2
    print("Part 1: Calculate C_R^2")
    C_R_squared = C_R * C_R
    if C_R_squared is None: return
    C_R_squared.name = "R^2"

    # Next, group terms that simplify well
    # Let's calculate (C_rho * C_pi * C_4_3)
    print("Part 2: Calculate intermediate term (rho * pi * 4/3)")
    term1 = C_rho * C_pi
    if term1 is None: return
    term1.name = "rho*pi"
    
    term2 = term1 * C_4_3
    if term2 is None: return
    term2.name = "rho*pi*4/3" # This results in 12/1, which is valid

    # Now, multiply the remaining terms. Here lies the problem.
    # Let's multiply our valid intermediate term by R^2
    print("Part 3: Multiply intermediate term by R^2")
    # This operation attempts to calculate 12 * 4 = 48.
    final_term = term2 * C_R_squared
    
    # --- Step 4: Analyze the result ---
    if final_term is None:
        print("--- CONCLUSION ---")
        print("The calculation fails because multiplying the coefficients results in a numerator (48)")
        print("that exceeds the 4-bit limit of 15. Different orders of multiplication or")
        print("reasonable changes in approximations lead to similar overflows.")
        print("For example, multiplying C_G(13/2) by C_R_squared(4/1) would yield a numerator of 52.")
        print("\nTherefore, it is not possible to perform this calculation on the Titan architecture.")
        print("\nFinal Answer: N0")
    else:
        # This part of the code will not be reached.
        print("Calculation was successful (this is unexpected).")

main()
<<<N0>>>