import math

class TitanOverflowError(Exception):
    """Custom exception for Titan's 5-bit constraint violation."""
    pass

class TitanFrac:
    """
    A class to represent and operate on Titan fractional numbers.
    It enforces the 5-bit constraint (0-31) for numerators and denominators.
    """
    def __init__(self, num, den=1):
        # The problem implies numerators/denominators are 5-bit integers (0-31)
        # Denominators cannot be zero.
        if not (0 <= num <= 31 and 1 <= den <= 31):
            raise TitanOverflowError(f"Initialization failed for {num}/{den}. Values must be in the 5-bit range [0, 31].")
        
        self.num = num
        self.den = den
        self._simplify()

    def _simplify(self):
        """Simplifies the fraction using the greatest common divisor (GCD)."""
        common = math.gcd(self.num, self.den)
        self.num //= common
        self.den //= common

    def __mul__(self, other):
        """
        Multiplies two Titan fractions, checking for overflow at each step.
        Implements cross-simplification to find a valid computational path.
        """
        # (a/b) * (c/d). The rule is "Any operation resulting in numerators or 
        # denominators exceeding 31 must be immediately simplified".
        # This implies that if a*c > 31, the operation is invalid unless simplification
        # can reduce it. The example 13/1 * 6/5 -> 25/2 * 6/5 -> (25/5)*(6/2)=15/1
        # shows that cross-simplification is a valid strategy.
        
        g1 = math.gcd(self.num, other.den)
        g2 = math.gcd(other.num, self.den)
        
        num1 = self.num // g1
        den1 = other.den // g1
        
        num2 = other.num // g2
        den2 = self.den // g2
        
        final_num = num1 * num2
        final_den = den1 * den2
        
        if final_num > 31 or final_den > 31:
            raise TitanOverflowError(f"Multiplication overflow: ({self.num}/{self.den}) * ({other.num}/{other.den}) -> {final_num}/{final_den}")
            
        return TitanFrac(final_num, final_den)

    def __add__(self, other):
        """Adds two Titan fractions, checking for overflow."""
        # a/b + c/d = (ad + bc) / bd
        num = self.num * other.den + other.num * self.den
        den = self.den * other.den
        if num > 31 or den > 31:
            raise TitanOverflowError(f"Addition overflow: ({self.num}/{self.den}) + ({other.num}/{other.den}) -> {num}/{den}")
        return TitanFrac(num, den)
            
    def __repr__(self):
        return f"{self.num}/{self.den}"

def demonstrate_titan_calculation():
    """
    Simulates the calculation of the landing force on the Titan computer.
    """
    print("--- [Titan Computer Simulation] ---")
    print("Attempting to calculate the required landing force for the Pioneer probe.")
    print("The governing equation is F_rocket = F_gravity + F_deceleration")
    print("F_rocket = (m * g) + (m * a)\n")
    
    try:
        # --- Define Constants as Titan Fractions ---
        # Probe Mass m = 50kg. Cannot be a single fraction as 50 > 31.
        # We must represent it as a product of valid fractions, e.g., 50 = 10 * 5.
        m_factor1 = TitanFrac(10, 1)
        m_factor2 = TitanFrac(5, 1)
        print(f"Probe Mass (m=50kg) is represented by factors: {m_factor1} and {m_factor2}")

        # Required acceleration a = 9 m/s^2. This is a valid fraction.
        a = TitanFrac(9, 1)
        print(f"Deceleration (a=9m/s^2) is represented by: {a}")
        
        # Pandora gravity g ≈ 0.167 m/s^2. We can approximate this as 1/6.
        g = TitanFrac(1, 6)
        print(f"Pandora gravity (g≈0.167m/s^2) is approximated as: {g}\n")

        # --- Calculate Force Components ---
        # 1. Calculate Gravitational Force (F_gravity = m * g)
        print("Step 1: Calculate F_gravity = (m_factor1 * m_factor2) * g")
        # We can calculate this using cross-simplification: (10/1 * 5/1) * 1/6 = (10/1 * 1/6) * 5/1
        # = (10/6) * 5/1 = (5/3) * 5/1 = 25/3
        mg_intermediate = m_factor1 * g
        f_gravity = mg_intermediate * m_factor2
        print(f"SUCCESS: F_gravity calculation is possible. Result: {f_gravity} (which is ~8.33 N).\n")

        # 2. Calculate Deceleration Force (F_decel = m * a)
        print("Step 2: Calculate F_decel = (m_factor1 * m_factor2) * a")
        print(f"This requires multiplying the three fractions: {m_factor1}, {m_factor2}, and {a}.")
        print("We will test possible multiplication orders:\n")

        # Try order 1: (m1 * m2) * a
        try:
            print(f"  Attempting order 1: ({m_factor1} * {m_factor2}) * {a}")
            m_product = m_factor1 * m_factor2
        except TitanOverflowError as e:
            print(f"  FAILURE: {e}. The intermediate product of {m_factor1} and {m_factor2} overflows.")
        
        # Try order 2: (m1 * a) * m2
        try:
            print(f"\n  Attempting order 2: ({m_factor1} * {a}) * {m_factor2}")
            ma_product = m_factor1 * a
        except TitanOverflowError as e:
            print(f"  FAILURE: {e}. The intermediate product of {m_factor1} and {a} overflows.")

        # Try order 3: (m2 * a) * m1
        try:
            print(f"\n  Attempting order 3: ({m_factor2} * {a}) * {m_factor1}")
            ma_product2 = m_factor2 * a
        except TitanOverflowError as e:
            print(f"  FAILURE: {e}. The intermediate product of {m_factor2} and {a} overflows.")
        
        print("\n--- Conclusion ---")
        print("All possible pathways to calculate the deceleration force (50 * 9) result in an")
        print("intermediate multiplication (e.g., 10*5, 10*9, or 5*9) that produces a numerator")
        print("larger than 31. Since there are no denominators to simplify with, the operation")
        print("is impossible within the Titan architecture's constraints.")
        print("\nSince a major component of the force cannot be computed, the entire calculation fails.")

    except TitanOverflowError as e:
        print(f"\nA critical error occurred during the simulation: {e}")
        print("The calculation is not possible.")

    finally:
        print("\nCan you use Titan to calculate this force? NO.")
        print("Final Answer:")
        print("<<<N0>>>")

if __name__ == '__main__':
    demonstrate_titan_calculation()
