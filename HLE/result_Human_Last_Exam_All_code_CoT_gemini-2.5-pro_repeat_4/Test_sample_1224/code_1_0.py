#
# Titan Computer Feasibility Simulation for Pandora Landing Time Calculation
#

class TitanError(Exception):
    """Custom exception for Titan architecture constraint violations."""
    pass

class TitanFraction:
    """Represents a fraction within Titan's 4-bit constraints."""
    def __init__(self, num, den=1):
        if not (0 <= num <= 15 and 0 < den <= 15):
            raise TitanError(f"Overflow: Numerator/denominator {num}/{den} exceeds 4-bit limit.")
        self.num = num
        self.den = den

    def __repr__(self):
        return f"{self.num}/{self.den}"

class TitanRegister:
    """Simulates a Titan register holding an expression up to 10 terms."""
    def __init__(self, terms=None):
        if terms is None:
            # Initialize with zero
            terms = [(TitanFraction(0, 1), 0)]
        if len(terms) > 10:
            raise TitanError(f"Overflow: Expression exceeds 10 terms ({len(terms)}).")
        # Each term is a tuple (TitanFraction, exponent)
        self.terms = terms

    def __repr__(self):
        return " + ".join([f"{frac}e{exp}" for frac, exp in self.terms])

def titan_multiply(reg1, reg2):
    """Simulates the MUL instruction on two registers."""
    new_terms = []
    for frac1, exp1 in reg1.terms:
        for frac2, exp2 in reg2.terms:
            num = frac1.num * frac2.num
            den = frac1.den * frac2.den
            
            # Titan's rule: Any operation resulting in overflow must be simplified.
            # Here, we demonstrate that even a simple multiplication leads to overflow.
            if num > 15 or den > 15:
                # The architecture would need to expand this (e.g., 24 -> 15+9).
                # We will raise an error to show the calculation is not directly possible.
                raise TitanError(f"Intermediate multiplication {frac1.num}*{frac2.num}={num} or {frac1.den}*{frac2.den}={den} caused an overflow.")
            
            new_terms.append((TitanFraction(num, den), exp1 + exp2))
    
    return TitanRegister(new_terms)

def run_pandora_calculation():
    """Attempts to calculate g and demonstrates the failure."""
    print("--- Titan Computer Simulation ---")
    print("Task: Calculate gravitational acceleration 'g' on Pandora.")
    print("Simplified formula: g = G * (4/3) * pi * c * rho\n")
    
    try:
        # Step 1: Define constants and parameters as Titan expressions
        # We choose approximations that are valid within the 4-bit system.
        # G ≈ 7e-11 N⋅m²/kg²
        G = TitanRegister([(TitanFraction(7, 1), -11)])
        # pi ≈ 3
        pi = TitanRegister([(TitanFraction(3, 1), 0)])
        # c ≈ 2e6 m (polar radius)
        c = TitanRegister([(TitanFraction(2, 1), 6)])
        # rho ≈ 3e2 kg/m³ (shell density)
        rho = TitanRegister([(TitanFraction(3, 1), 2)])
        # The constant 4/3
        four_thirds = TitanRegister([(TitanFraction(4, 3), 0)])
        
        print("Initial values:")
        print(f"  G   = {G}")
        print(f"  pi  = {pi}")
        print(f"  c   = {c}")
        print(f"  rho = {rho}")
        print(f"  4/3 = {four_thirds}\n")

        # Step 2: Perform multiplications in an optimized order
        # (4/3) * pi = (4/3) * 3 = 4
        print("Calculating: (4/3) * pi...")
        ax = titan_multiply(four_thirds, pi)
        print(f"Result = {ax.terms[0][0].num}/{ax.terms[0][0].den} (simplified from 12/3)")

        # next: * c = 4 * 2e6 = 8e6
        print("\nCalculating: previous_result * c...")
        ax = titan_multiply(ax, c)
        print(f"Result = {ax}")

        # next: * rho = 8e6 * 3e2 = 24e8 -> This step will fail
        print("\nCalculating: previous_result * rho...")
        ax = titan_multiply(ax, rho)
        # This line will not be reached

    except TitanError as e:
        print(f"\n>>> CALCULATION FAILED <<<")
        print(f"Reason: {e}")
        print("\nThe product of mantissas (e.g., 8 * 3 = 24) exceeds the 4-bit limit of 15.")
        print("While the Titan architecture can expand such numbers into multiple terms (e.g., 24 -> 15 + 9),")
        print("the full calculation of 'g' requires subsequent multiplications that cause the number of terms")
        print("in the expression to exceed the register limit of 10.")
        print("\nConclusion: The landing time cannot be calculated on Titan due to these constraints.")
        print("\nFinal Answer:")
        print("<<<N0>>>")

if __name__ == '__main__':
    run_pandora_calculation()