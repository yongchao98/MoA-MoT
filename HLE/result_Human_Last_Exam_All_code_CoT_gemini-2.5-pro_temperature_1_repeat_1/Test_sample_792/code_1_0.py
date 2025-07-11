import math

def gcd(a, b):
    """Helper function to find the greatest common divisor."""
    return math.gcd(a, b)

class TitanFraction:
    """
    A class to represent and compute with Titan's 5-bit fractional numbers.
    It automatically handles simplification and validation against 5-bit constraints.
    """
    def __init__(self, num, den=1):
        if den == 0:
            raise ValueError("Denominator cannot be zero.")
        
        # Simplify the fraction upon creation
        common = gcd(num, den)
        self.num = num // common
        self.den = den // common

        # Check if the simplified fraction is within 5-bit limits
        if not (0 <= self.num <= 31 and 1 <= self.den <= 31):
            raise ValueError(
                f"Value {num}/{den} (simplified to {self.num}/{self.den}) is out of the 5-bit integer range (0-31)."
            )

    def __mul__(self, other):
        """Multiplies two TitanFractions, applying simplification rules."""
        print(f"Multiplying: ({self.num}/{self.den}) * ({other.num}/{other.den})")
        new_num = self.num * other.num
        new_den = self.den * other.den
        print(f"  - Raw intermediate result: {new_num}/{new_den}")
        # The constructor will handle simplification and validation
        result = TitanFraction(new_num, new_den)
        print(f"  - Valid simplified result: {result}")
        return result

    def __truediv__(self, other):
        """Divides two TitanFractions by multiplying by the reciprocal."""
        print(f"Dividing: ({self.num}/{self.den}) / ({other.num}/{other.den})")
        if other.num == 0:
            raise ZeroDivisionError("Titan fraction division by zero")
        
        # The constructor handles the reciprocal's creation and validation
        reciprocal = TitanFraction(other.den, other.num)
        print(f"  - Rewriting as multiplication: ({self.num}/{self.den}) * ({reciprocal})")
        return self * reciprocal

    def __str__(self):
        return f"{self.num}/{self.den}"

def calculate_force_on_titan():
    """
    Simulates the entire calculation on the Titan computer.
    """
    try:
        print("--- Titan Superconducting Computer Booting Up ---")
        print("Objective: Calculate the force required to hit the lion.")
        print("-" * 50)

        # Step 1: Define constants and approximations as Titan Fractions.
        print("Step 1: Loading constants and approximations...")
        h = TitanFraction(10, 1)          # Height = 10m
        d_target = TitanFraction(20, 1)   # Target distance = 20m
        r = TitanFraction(1, 2)           # Rock radius = 0.5cm
        four_thirds = TitanFraction(4, 3)

        # Approximations chosen to be computationally viable on Titan.
        g_approx = TitanFraction(10, 1)   # g ≈ 10 m/s^2 (from 9.8)
        pi_approx = TitanFraction(3, 1)   # π ≈ 3 (from 3.14159...)
        rho_approx = TitanFraction(9, 10) # Rock density = 0.9 kg/cm^3

        print(f"h = {h}, d = {d_target}, r = {r}, 4/3 = {four_thirds}")
        print(f"Approximations: g ≈ {g_approx}, π ≈ {pi_approx}, ρ ≈ {rho_approx}")
        print("-" * 50)

        # Step 2: Calculate the mass (m) of the rock. m = ρ * (4/3) * π * r³
        print("Step 2: Calculating rock mass (m)...")
        r_squared = r * r
        r_cubed = r_squared * r
        
        m_intermediate1 = rho_approx * four_thirds
        m_intermediate2 = m_intermediate1 * pi_approx
        m = m_intermediate2 * r_cubed
        print(f"Final calculated mass: m = {m} kg")
        print("-" * 50)

        # Step 3: Calculate the required force (F). F = (d/h) * m * g
        print("Step 3: Calculating required force (F)...")
        print("Re-ordering to F = (d / h) * m * g for numerical stability.")
        
        d_over_h = d_target / h
        f_intermediate = d_over_h * m
        F = f_intermediate * g_approx
        print("-" * 50)

        # Step 4: Final Result
        print("Final Result from Titan Calculation:")
        print(f"The equation for the force is F = (d * m * g) / h.")
        print(f"The required force is: {F.num} / {F.den} N")

    except (ValueError, ZeroDivisionError) as e:
        print(f"\n--- COMPUTATION FAILED ---")
        print(f"An error occurred during calculation: {e}")

# Run the simulation
calculate_force_on_titan()