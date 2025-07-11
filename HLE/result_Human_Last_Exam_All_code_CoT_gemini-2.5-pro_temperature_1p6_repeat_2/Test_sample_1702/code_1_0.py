import math

class WuxingFrac:
    """
    A Python class to simulate the Wuxing 'frac' data type.
    It holds a value as (n/d) * 10^e.
    """
    def __init__(self, n, d, e):
        if d == 0:
            raise ZeroDivisionError("Denominator cannot be zero.")
        
        # Automatically simplify the fraction n/d using the greatest common divisor
        common = math.gcd(n, d)
        self.n = n // common
        self.d = d // common
        self.e = e

    def __mul__(self, other):
        # Handle multiplication with integers (like the '2' in the formula)
        if isinstance(other, int):
            other = WuxingFrac(other, 1, 0)
        
        new_n = self.n * other.n
        new_d = self.d * other.d
        new_e = self.e + other.e
        return WuxingFrac(new_n, new_d, new_e)
    
    # Required for right-side multiplication (e.g., 2 * frac_object)
    def __rmul__(self, other):
        return self.__mul__(other)

    def __truediv__(self, other):
        # Handle division by integers (like the final division by 1000)
        if isinstance(other, int):
            other = WuxingFrac(other, 1, 0)

        new_n = self.n * other.d
        new_d = self.d * other.n
        new_e = self.e - other.e
        return WuxingFrac(new_n, new_d, new_e)

    def __pow__(self, power):
        if not isinstance(power, int) or power < 0:
            raise ValueError("Power must be a non-negative integer.")
        result = WuxingFrac(1, 1, 0)
        for _ in range(power):
            result = result * self
        return result

    def get_value(self):
        """Returns the standard decimal value of the fraction."""
        return (self.n / self.d) * (10 ** self.e)

def main():
    """
    Calculates Pandora's temperature based on the Rayleigh-Jeans Law,
    simulating the Wuxing architecture's fractional arithmetic.
    """
    # Define constants and given values as WuxingFrac objects.
    # Note: Using k ≈ 1.4e-23 to respect the `signed char` numerator constraint (n=14 is valid).
    B_lambda = WuxingFrac(99, 10, 16)   # Spectral Radiance B = 9.9e16
    lambda_ = WuxingFrac(5, 1, -7)      # Wavelength λ = 5e-7 m
    c = WuxingFrac(3, 1, 8)             # Speed of light c = 3e8 m/s
    k = WuxingFrac(14, 10, -23)         # Boltzmann constant k ≈ 1.4e-23 J/K
    
    # Calculate T = (B * λ^4) / (2 * c * k)
    lambda_4 = lambda_ ** 4
    numerator = B_lambda * lambda_4
    denominator = 2 * c * k
    temperature_K = numerator / denominator

    # The result needs to be in thousands of Kelvin.
    temperature_kilo_K = temperature_K / 1000

    # Get the decimal value and round it to the nearest whole number.
    final_value = temperature_kilo_K.get_value()
    rounded_temp_in_kilo_K = round(final_value)
    
    # As requested, print the final equation with all its numbers.
    print("Final Equation: T = (9.9e16 * (5e-7)^4) / (2 * 3e8 * 1.4e-23)")
    print(rounded_temp_in_kilo_K)

if __name__ == "__main__":
    main()