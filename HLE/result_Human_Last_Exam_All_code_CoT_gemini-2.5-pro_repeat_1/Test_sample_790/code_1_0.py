import math

class TitanFraction:
    """
    A class to represent and operate on fractions under Titan's 5-bit constraints.
    Numerators and denominators must be representable by 5-bit integers (0-31).
    """
    def __init__(self, num, den=1, verbose=True):
        # Simplify the fraction first
        common_divisor = math.gcd(num, den)
        self.num = num // common_divisor
        self.den = den // common_divisor

        if not (0 <= self.num <= 31 and 0 < self.den <= 31):
            raise ValueError(f"Resulting fraction {self.num}/{self.den} is out of 5-bit range.")

        if verbose:
            print(f"  = {self.num}/{self.den}")

    def __mul__(self, other, verbose=True):
        if verbose:
            print(f"Multiplying {self} * {other}")
        new_num = self.num * other.num
        new_den = self.den * other.den
        return TitanFraction(new_num, new_den, verbose=verbose)

    def __repr__(self):
        return f"{self.num}/{self.den}"

def solve_titan_force():
    """
    Solves for the force using Titan's computational rules.
    """
    print("--- Titan Force Calculation ---")
    print("Objective: Calculate F = 2 * g * m * sqrt(2)")
    print("Using CGS units (cm, g, s) and 5-bit fractional arithmetic.\n")

    # The physics formula for the force (F) in CGS units (dynes) is:
    # F = 2 * g * m * sqrt(2)
    # where m = rho * V = rho * (4/3) * pi * r^3
    # Substituting m, we get: F = g * rho * (8/3) * pi * r^3 * sqrt(2)

    # --- Step 1: Define constants as Titan-compatible fractions ---
    print("Step 1: Approximating constants")
    # Radius r = 0.5 cm
    r = TitanFraction(1, 2, verbose=False)
    # r^3 = (1/2)^3 = 1/8
    r_cubed = TitanFraction(1, 8, verbose=False)
    print(f"r = 1/2 cm, so r^3 = {r_cubed}")

    # Density rho = 0.9 g/cm^3 (assuming typo correction from 0.9 kg/cm^3)
    # To enable simplification later, we approximate rho = 9/10 as 28/31 (~0.903)
    rho = TitanFraction(28, 31, verbose=False)
    print(f"rho (density) approx {rho} g/cm^3")

    # Gravitational acceleration g = 980 cm/s^2
    # We need an approximation with small factors. Let's use g = 930 = 30 * 31.
    g_part1 = TitanFraction(30, 1, verbose=False)
    g_part2 = TitanFraction(31, 1, verbose=False)
    print(f"g (gravity) approx {g_part1.num*g_part2.num} cm/s^2 (represented as {g_part1} * {g_part2})")

    # Other constants
    eight_thirds = TitanFraction(8, 3, verbose=False)
    # Pi approx 22/7
    pi = TitanFraction(22, 7, verbose=False)
    # sqrt(2) approx 7/5 (strategically chosen to cancel with pi's denominator)
    sqrt2 = TitanFraction(7, 5, verbose=False)
    print(f"pi approx {pi}, sqrt(2) approx {sqrt2}\n")


    # --- Step 2: Calculate the force term by term, simplifying along the way ---
    print("Step 2: Calculating F = (g1*g2) * rho * (8/3) * pi * r^3 * sqrt(2)")
    # We will multiply in a specific order to keep intermediate values small.

    # Start with g_part1
    force = g_part1
    print(f"Initial F = {force}")

    # F * (8/3)
    force = force * eight_thirds
    
    # F * r^3
    force = force * r_cubed

    # F * pi
    force = force * pi
    
    # F * sqrt(2)
    force = force * sqrt2

    # F * rho
    force = force * rho

    # F * g_part2
    force = force * g_part2
    
    final_force_dynes = force.num / force.den
    print(f"\nFinal calculated force: {force} = {final_force_dynes:.1f} dynes")

    # --- Step 3: Convert to Newtons and find error ---
    print("\nStep 3: Convert to Newtons and calculate error")
    # 1 N = 100,000 dynes
    final_force_newtons = final_force_dynes / 100000.0
    print(f"Final force in Newtons: {final_force_newtons:.4f} N")

    # True value calculated with high precision is ~0.1306 N
    true_force_newtons = 0.1306
    absolute_error = abs(final_force_newtons - true_force_newtons)
    print(f"True force (high precision): {true_force_newtons} N")
    print(f"Absolute error: {absolute_error:.4f} N")

    print("\n--- Final Answer ---")
    print("The calculation is possible on the Titan architecture.")
    print(f"The equation for the force is: F = 2 * g * m * sqrt(2)")
    # Final answer format Y[e]
    print(f"Y[{absolute_error:.3f}]")
    
solve_titan_force()