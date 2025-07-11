import math

class TitanComputer:
    """
    A class to simulate the Titan computer's fractional arithmetic.
    It enforces the 5-bit constraint (numerators and denominators <= 31).
    """

    @staticmethod
    def create_fraction(n, d=1, check_overflow=True):
        """Creates a fraction if it meets Titan's constraints."""
        if check_overflow and not (0 <= n <= 31 and 1 <= d <= 31):
            # This is a fundamental constraint violation.
            # We return None to signify that the number cannot be represented.
            return None
        
        common = math.gcd(n, d)
        return (n // common, d // common)

    def multiply(self, f1, f2):
        """Multiplies two fractions, simplifying before to prevent overflow."""
        if f1 is None or f2 is None: return None
        n1, d1 = f1
        n2, d2 = f2
        
        # Pre-simplify to avoid intermediate overflow, as shown in the example
        g1 = math.gcd(n1, d2)
        g2 = math.gcd(n2, d1)
        
        num = (n1 // g1) * (n2 // g2)
        den = (d1 // g2) * (d2 // g1)
        
        return self.create_fraction(num, den)

    def add(self, f1, f2):
        """Adds two fractions."""
        if f1 is None or f2 is None: return None
        n1, d1 = f1
        n2, d2 = f2
        
        num = n1 * d2 + n2 * d1
        den = d1 * d2
        
        return self.create_fraction(num, den)

def solve():
    """
    Main function to attempt the gravitational force calculation.
    """
    titan = TitanComputer()

    # Step 1: Define physical and mathematical constants as Titan fractions.
    # We will analyze each constant's feasibility.
    
    # Probe mass: m_p = 30 kg. 30 is a valid 5-bit integer.
    m_p = titan.create_fraction(30, 1)
    
    # Pi: approx 22/7. Both 22 and 7 are valid 5-bit integers.
    pi = titan.create_fraction(22, 7)
    
    # 4/3: a constant from the volume formula for a sphere.
    four_thirds = titan.create_fraction(4, 3)

    # Now, let's analyze the planetary data. We use SI units for consistency.
    # Shell radius: R_p = 1000 km = 1,000,000 m
    # Shell density: ρ_s = 0.3 metric tons/m^3 = 300 kg/m^3
    # Gravitational constant: G ≈ 6.674e-11 N m^2/kg^2
    
    # The values 1,000,000 (for R_p) and 300 (for ρ_s) are vastly larger
    # than the maximum 5-bit integer (31). They cannot be represented.
    R_p = titan.create_fraction(1000000, 1) # This will fail
    rho_s = titan.create_fraction(300, 1) # This will fail
    
    # The constant G is also impossible to represent due to its small magnitude (10^-11).
    
    # The problem cannot be solved by direct computation with the given values.
    # An alternative is to simplify the physics equation algebraically to see if
    # these large numbers cancel out or form representable ratios.
    
    # Simplified Force Equation: F ≈ G * (4/3)π * ρ_s * m_p * R_p
    # This simplification requires two assumptions that are forced by Titan's limits:
    # 1. Altitude `h` is negligible compared to `R_p` (h/R_p = 500/1000000 = 1/2000, denominator too large)
    # 2. The core's extra mass contribution is negligible (3*(r_c/R_p)^3 is a small term that would cause overflow)
    
    # This simplified equation still contains the product: G * ρ_s * R_p
    # Let's calculate its value: 6.674e-11 * 300 * 1,000,000 ≈ 0.02
    # To be represented on Titan, 0.02 would be 2/100 = 1/50.
    # The denominator 50 is larger than 31, so this value also cannot be represented.

    # Conclusion: The fundamental constants of the problem (planet size, density, G)
    # have magnitudes that fall outside the representational capacity of the Titan computer.
    # The constraints are too restrictive to allow for a calculation.
    if R_p is None or rho_s is None:
        print("N0")
    else:
        # This part of the code is unreachable because the constants are not representable.
        # If they were, the calculation would proceed here.
        pass

solve()
<<<N0>>>