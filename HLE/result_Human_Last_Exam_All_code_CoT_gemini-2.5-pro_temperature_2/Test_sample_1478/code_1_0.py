import sys

# This program simulates the calculation on the Titan 6-bit architecture.
# We will demonstrate that the calculation is not possible due to the constraints.

# Titan computer's maximum integer value for numerators and denominators.
# A 6-bit unsigned integer can represent values from 0 to 2^6 - 1 = 63.
MAX_VAL = 63

class TitanTerm:
    """Represents a single term (fraction * 10^exp) in a Titan register."""
    def __init__(self, num, den=1, exp=0):
        if abs(num) > MAX_VAL or abs(den) > MAX_VAL:
            # This class helps us track when a value goes out of bounds.
            raise ValueError(f"Overflow Error: {num}/{den} exceeds 6-bit limit.")
        self.num = num
        self.den = den
        self.exp = exp

    def __str__(self):
        # We only print the fraction part for clarity, as exponents are handled separately.
        return f"{self.num}/{self.den}"

def print_step(step_num, instruction, result, comment):
    """Helper function to format the output."""
    print(f"Step {step_num}: {instruction:<15} | AX = {result:<25} | {comment}")

def check_overflow(num, den=1):
    """Checks if a value would overflow."""
    return abs(num) > MAX_VAL or abs(den) > MAX_VAL

def main():
    print("--- Starting Titan Feasibility Analysis ---")
    print(f"Constraint: All numerators/denominators must be <= {MAX_VAL}\n")

    # Step 1: Define physical quantities as Titan fractions.
    # We will focus on the mantissas, as exponents just add up.
    print("1. Define constants and initial values:")
    rho_m = TitanTerm(12, 1)      # Density of Pandora (1200 -> 12e2)
    four_thirds = TitanTerm(4, 3) # From volume formula
    pi_approx = TitanTerm(22, 7)  # Approximation of Pi
    r_cubed_m = TitanTerm(8, 1)   # Radius of Pandora cubed (2^3=8)
    m_probe = TitanTerm(50, 1)    # Probe Mass
    G_m = TitanTerm(20, 3)        # Gravitational Constant (6.66.. -> 20/3)
    d_squared_m = TitanTerm(1,1)  # Distance from black hole squared (1km=1000m -> 1e3, mantissa is 1)
    
    print(f"   ρ_m = {rho_m}, 4/3 = {four_thirds}, π_m = {pi_approx}, r³_m = {r_cubed_m}")
    print("-" * 50)
    print("2. Attempt to calculate Mass mantissa M_m = ρ_m * 4/3 * π_m * r³_m")
    
    # --- Simulation of Mass Calculation ---
    
    # MOV AX, rho_m
    ax_terms = [rho_m]
    print_step(1, "MOV AX, 12/1", "12/1", "Load density mantissa.")
    
    # MUL AX, 4/3
    term1 = ax_terms[0]
    term2 = four_thirds
    new_num = term1.num * term2.num # 12 * 4 = 48
    new_den = term1.den * term2.den # 1 * 3 = 3
    # Simplify the fraction by finding the greatest common divisor (GCD). Here, GCD(48,3)=3
    simplified_num = new_num // 3
    simplified_den = new_den // 3
    ax_terms = [TitanTerm(simplified_num, simplified_den)] # Result is 16/1
    print_step(2, "MUL AX, 4/3", "16/1", "OK. Result is within limits.")

    # MUL AX, 22/7
    term1 = ax_terms[0]
    term2 = pi_approx
    # This direct multiplication would result in 16*22=352, causing an overflow.
    # We must use the reduction rule: 16 * (3 + 1/7) = 16*3 + 16*7
    op1_num = term1.num * 3
    op1_den = term1.den
    op2_num = term1.num * 1
    op2_den = term1.den * 7

    if check_overflow(op1_num, op1_den) or check_overflow(op2_num, op2_den):
        # This check is technically not needed here as 48 and 16/7 are valid.
        pass
    
    ax_terms = [TitanTerm(op1_num, op1_den), TitanTerm(op2_num, op2_den)]
    print_step(3, "MUL AX, 22/7", "[48/1, 16/7]", "Overflow avoided via reduction. AX holds an expression.")

    # MUL AX, 8/1
    # Now we must multiply the expression [48/1, 16/7] by 8/1.
    # This becomes (48/1 * 8/1) + (16/7 * 8/1)
    term_a = ax_terms[0]
    term_b = ax_terms[1]
    r_term = r_cubed_m

    final_term1_num = term_a.num * r_term.num # 48 * 8 = 384
    final_term2_num = term_b.num * r_term.num # 16 * 8 = 128
    
    print_step(4, "MUL AX, 8/1", f"[{final_term1_num}/1, {final_term2_num}/7]", "Attempt to multiply expression by r³_m.")

    print("\n--- ANALYSIS OF STEP 4 ---")
    if check_overflow(final_term1_num) or check_overflow(final_term2_num):
        print(f"FAILURE: The operation resulted in the expression [{final_term1_num}/1, {final_term2_num}/7].")
        print(f"The numerator of the first term ({final_term1_num}) exceeds {MAX_VAL}.")
        print(f"The numerator of the second term ({final_term2_num}) exceeds {MAX_VAL}.")
        print("\nCONCLUSION: The mass M cannot be calculated on the Titan architecture.")
        print("The calculation is impossible with the given constraints.")
        print("Therefore, the force cannot be determined.")

if __name__ == "__main__":
    main()
