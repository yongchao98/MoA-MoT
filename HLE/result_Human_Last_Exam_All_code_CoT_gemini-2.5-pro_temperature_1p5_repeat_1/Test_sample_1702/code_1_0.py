import math

class Frac:
    """
    A class to simulate the Wuxing 'frac' type.
    It represents a number as (n/d) * 10^e.
    """
    def __init__(self, n, d, e):
        if d == 0:
            raise ValueError("Denominator cannot be zero.")
        self.n = int(n)
        self.d = int(d)
        self.e = int(e)
        self._simplify()

    def _simplify(self):
        """Simplifies the fraction n/d."""
        if self.n == 0:
            self.d = 1
            return
        common = math.gcd(self.n, self.d)
        self.n //= common
        self.d //= common
        if self.d < 0:
            self.n = -self.n
            self.d = -self.d

    def __mul__(self, other):
        """Multiplies two Frac numbers."""
        n = self.n * other.n
        d = self.d * other.d
        e = self.e + other.e
        return Frac(n, d, e)

    def __truediv__(self, other):
        """Divides two Frac numbers."""
        n = self.n * other.d
        d = self.d * other.n
        e = self.e - other.e
        return Frac(n, d, e)

    def to_float(self):
        """Converts the Frac number to a standard float."""
        return (self.n / self.d) * (10 ** self.e)

    def __str__(self):
        """String representation in the format n/d*10^e."""
        return f"{self.n}/{self.d}e{self.e}"

def solve_pandora_temperature():
    """
    Calculates the temperature of Pandora using the Rayleigh-Jeans approximation
    and the custom Frac type.
    """
    # Given values and constants as Frac objects
    # Spectral Radiance Bλ = 9.9e16 W/m^2·sr·m
    B_lambda = Frac(99, 10, 16)
    # Wavelength λ = 500 nm = 5e-7 m
    lambda_ = Frac(5, 1, -7)
    # Constant 2
    two = Frac(2, 1, 0)
    # Speed of light c ≈ 3e8 m/s
    c = Frac(3, 1, 8)
    # Boltzmann constant kB ≈ 1.38e-23 J/K. Approximated as 11/8 = 1.375.
    kB = Frac(11, 8, -23)

    # Temperature Formula: T = (Bλ * λ⁴) / (2 * c * kB)

    # Calculate λ⁴
    lambda_2 = lambda_ * lambda_
    lambda_4 = lambda_2 * lambda_2

    # Calculate the full numerator and denominator Frac objects
    numerator_full = B_lambda * lambda_4
    denominator_full = two * c * kB

    # As derived in the thinking process, the simplified fraction forms are:
    # Numerator: 99/16e-9
    # Denominator: 33/4e-15
    
    # Calculate the temperature and the final rounded answer
    temperature_K = (numerator_full / denominator_full).to_float()
    final_answer = round(temperature_K / 1000)
    
    # Print the final equation with the simplified fractional numbers
    print("Final Equation:")
    print("T = (99/16e-9) / (33/4e-15)")
    print(f"Calculated Temperature: {temperature_K} K")
    print(f"Temperature in thousand Kelvin (rounded): {final_answer}")


solve_pandora_temperature()
<<<750>>>