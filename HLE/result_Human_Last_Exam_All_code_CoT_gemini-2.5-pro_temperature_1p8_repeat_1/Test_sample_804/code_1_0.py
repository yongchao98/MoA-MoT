import math

class TitanComputerError(Exception):
    """Custom exception for Titan computation errors."""
    pass

class Fraction:
    """A class to represent fractions under Titan's 5-bit constraints."""
    def __init__(self, numerator, denominator=1):
        if not (0 <= numerator <= 31 and 1 <= denominator <= 31):
            raise TitanComputerError(f"Initialization failed: Numerator ({numerator}) or denominator ({denominator}) out of 5-bit range (0-31).")
        self.n = numerator
        self.d = denominator

    def __str__(self):
        return f"{self.n}/{self.d}"

    def value(self):
        return self.n / self.d

def titan_simplify(n, d):
    """Simplifies a fraction n/d by dividing by their greatest common divisor."""
    if d == 0:
        raise TitanComputerError("Division by zero.")
    common_divisor = math.gcd(n, d)
    return n // common_divisor, d // common_divisor

def titan_multiply(f1, f2):
    """Multiplies two Fractions, checking for overflow at each step."""
    # Cross-simplify first to minimize intermediate values
    n1, d2 = titan_simplify(f1.n, f2.d)
    n2, d1 = titan_simplify(f2.n, f1.d)
    
    # Check for overflow before final multiplication
    if n1 * n2 > 31 or d1 * d2 > 31:
        raise TitanComputerError(
            f"Overflow during multiplication of {f1} * {f2}. "
            f"Intermediate result would be {n1*n2}/{d1*d2}, which exceeds the 5-bit limit."
        )
    return Fraction(n1 * n2, d1 * d2)

def titan_add(f1, f2):
    """Adds two Fractions, checking for overflow."""
    # n1/d1 + n2/d2 = (n1*d2 + n2*d1) / (d1*d2)
    n_part1 = f1.n * f2.d
    n_part2 = f2.n * f1.d
    d_final = f1.d * f2.d
    
    if n_part1 > 31 or n_part2 > 31 or d_final > 31:
         raise TitanComputerError(f"Overflow during addition of {f1} + {f2}. Intermediate products are too large.")
         
    n_final = n_part1 + n_part2
    if n_final > 31:
        raise TitanComputerError(f"Overflow during addition of {f1} + {f2}. Final numerator {n_final} is too large.")

    n_simple, d_simple = titan_simplify(n_final, d_final)

    if n_simple > 31 or d_simple > 31:
        raise TitanComputerError("Overflow in addition even after simplification.")

    return Fraction(n_simple, d_simple)

def main():
    """
    Attempts to calculate Pandora's gravity and determines if it's feasible on Titan.
    """
    print("Step 1: Define physical quantities as ratios and fractions.")
    # Physical values
    # r_core = 50 km, r_total = 1000 km -> ratio = 50/1000 = 1/20
    # The simplest fractional representation for this ratio is r_c_frac being 1/1 and r_t_frac being 20/1.
    # However, to maintain the ratio as a fraction itself, we use r_ratio = 1/20.
    try:
        r_ratio = Fraction(1, 20)
        print(f"Ratio of radii (r_core/r_total): {r_ratio}")
    except TitanComputerError as e:
        print(f"Error: {e}")
        print("\nConclusion: The problem is not solvable.")
        print("<<<N0>>>")
        return

    print("\nStep 2: Formulate the mass calculation using these fractions.")
    # The mass of the planet is proportional to: M ∝ (r_total)³ * [ (r_core/r_total)³ * (ρ_core/ρ_shell - 1) + 1 ] * ρ_shell
    # This calculation requires cubing the radius ratio, i.e., r_ratio * r_ratio * r_ratio.
    print(f"A key step is to calculate the cube of the radius ratio: ({r_ratio})^3.")

    print("\nStep 3: Attempt the calculation with Titan's arithmetic rules.")
    try:
        print(f"Calculating ({r_ratio}) * ({r_ratio})...")
        r_ratio_sq = titan_multiply(r_ratio, r_ratio)
        # This step will fail because 1*1=1 (ok), but 20*20=400, which is > 31.
        print(f"Result of squaring: {r_ratio_sq}")
        r_ratio_cubed = titan_multiply(r_ratio_sq, r_ratio)
        print(f"Result of cubing: {r_ratio_cubed}")
    except TitanComputerError as e:
        print(f"Calculation failed.")
        print(f"Reason: {e}")
        print("\nThe gravitational force cannot be computed on the Titan architecture because fundamental operations,")
        print("like squaring the ratio of the planet's radii, result in numbers that exceed the 5-bit limit.")
        print("The architecture is unable to handle the scale of the numbers involved, even when expressed as ratios.")
        print("\nFinal symbolic equation (uncalculable): F = G * m * (4/3) * pi * (r_total^3 * rho_shell + r_core^3 * (rho_core - rho_shell)) / r^2")
        print("With r_core/r_total = 1/20, the term (1/20)^2 results in 1/400, which is not representable.")
        print("\n<<<N0>>>")

if __name__ == "__main__":
    main()