import math

class TitanError(ValueError):
    """Custom error for Titan's 4-bit constraint violation."""
    pass

class TitanNumber:
    """
    A class to represent a number within Titan's 4-bit fractional system.
    It enforces the constraint that numerators and denominators must be <= 15.
    """
    def __init__(self, num, den=1, silent=False):
        # The check must happen for any operation.
        if not (0 <= num <= 15 and 1 <= den <= 15):
            if not silent:
                 raise TitanError(f"Constraint Violation: The fraction {num}/{den} is invalid. Numerators and denominators must be between 0 and 15.")
        self.num = num
        self.den = den

    def __str__(self):
        if self.den == 1:
            return f"{self.num}"
        return f"{self.num}/{self.den}"

def main():
    """
    Analyzes the feasibility of the Pandora landing time problem on the Titan architecture.
    """
    print("--- Titan Computer Feasibility Analysis for Pandora Landing Time ---")
    print("\nObjective: Calculate landing time t = sqrt(2*d/a). This requires calculating gravitational acceleration 'a'.")
    print("Formula for gravity: a = G * M / R^2. The main task is to compute Pandora's mass, M.")

    print("\nStep 1: Simplify the problem and choose 4-bit approximations.")
    print("We approximate Pandora as a simple sphere and its mass as the mass of its shell (M ≈ M_shell), since the core's contribution is much smaller.")
    print("M_shell = ρ_shell * Volume = ρ_shell * (4/3) * π * R³")

    # Use simple integer approximations for mantissas to test the core logic
    # ρ_shell = 300 kg/m^3 -> mantissa 3
    # R ≈ 2000 km -> mantissa 2
    # π ≈ 3
    rho = TitanNumber(3)
    four_thirds = TitanNumber(4, 3)
    pi = TitanNumber(3)
    radius_mantissa = TitanNumber(2)
    radius_cubed_mantissa = TitanNumber(8) # 2*2*2 = 8

    print("\nThe equation for the mantissa of the mass is: M_mantissa = ρ * (4/3) * π * R³")
    print("Using our 4-bit integer approximations, the final equation is:")
    print(f"M_mantissa = {rho} * {four_thirds} * {pi} * {radius_cubed_mantissa}\n")

    print("Step 2: Simulate the calculation on Titan, step-by-step.")
    try:
        # MOV AX, 3/1 (rho)
        # MUL AX, 4/3
        print("Instruction 1: MUL 3/1, 4/3")
        res1_num = rho.num * four_thirds.num     # 3 * 4 = 12
        res1_den = rho.den * four_thirds.den     # 1 * 3 = 3
        print(f"  - Intermediate product: {res1_num}/{res1_den}. This is valid.")
        # Simplify the fraction. 12/3 simplifies to 4/1.
        current_result = TitanNumber(4, 1)
        print(f"  - Result in register after simplification: {current_result} (Valid)\n")

        # MUL AX, 3/1 (pi)
        print("Instruction 2: MUL 4/1, 3/1")
        res2_num = current_result.num * pi.num   # 4 * 3 = 12
        res2_den = current_result.den * pi.den   # 1 * 1 = 1
        print(f"  - Intermediate product: {res2_num}/{res2_den}. This is valid.")
        current_result = TitanNumber(res2_num, res2_den)
        print(f"  - Result in register: {current_result} (Valid)\n")

        # MUL AX, 8/1 (R^3)
        print("Instruction 3: MUL 12/1, 8/1")
        res3_num = current_result.num * radius_cubed_mantissa.num # 12 * 8 = 96
        res3_den = current_result.den * radius_cubed_mantissa.den # 1 * 1 = 1
        print(f"  - Intermediate product: {res3_num}/{res3_den}.")
        print("\n  - Checking constraint for the result...")
        # This final step will fail because the numerator 96 is > 15
        final_result = TitanNumber(res3_num, res3_den)

    except TitanError as e:
        print(f"  - {e}\n")

    print("Step 3: Conclusion")
    print("The calculation fails because the multiplication of two valid 4-bit numbers (12 and 8) produces a result (96) whose numerator exceeds the 4-bit limit of 15.")
    print("This issue is fundamental. The problem's physical constants are too large to be handled by the Titan architecture's strict multiplication rules.")
    print("Workarounds like expression expansion or aggressive rounding (RDX) lead to either similarly invalid results or a catastrophic loss of precision, rendering the final answer meaningless.")
    print("\nTherefore, it is not feasible to use the Titan computer for this task.")

if __name__ == "__main__":
    main()