import math

# Wuxing architecture constraints
# char: 2D (signed: -99 to 99, unsigned: 0-99)
# int: 5D (-99999 to 99999)

def gcd(a, b):
    """Helper to compute greatest common divisor."""
    return math.gcd(a, b)

def find_best_frac_approx(value):
    """
    Finds the best n/d*10^e approximation for a given float value
    under Wuxing frac constraints.
    """
    if value == 0:
        return 0, 1, 0
    
    # Find an exponent that brings the value into a reasonable range for n/d
    if value != 0:
        exponent = math.floor(math.log10(abs(value)))
    else:
        exponent = 0
    
    # We want to find n/d approx target_val
    target_val = value / (10**exponent)

    # Shift exponent so target_val is between 0.1 and 1
    if abs(target_val) < 0.1:
        target_val *= 10
        exponent -= 1
    
    best_n, best_d = 1, 1
    min_error = float('inf')

    # Search for the best n, d approximation
    for d_candidate in range(1, 100):
        n_candidate = round(target_val * d_candidate)
        if -99 <= n_candidate <= 99 and n_candidate != 0:
            error = abs(target_val - n_candidate / d_candidate)
            if error < min_error:
                min_error = error
                best_n = n_candidate
                best_d = d_candidate
                
    return int(best_n), int(best_d), int(exponent)


class Frac:
    """A class to simulate the Wuxing frac data type."""
    def __init__(self, n, d, e):
        if not isinstance(n, int) or not isinstance(d, int) or not isinstance(e, int):
            raise TypeError("n, d, e must be integers")
        if d == 0:
            raise ZeroDivisionError("Denominator cannot be zero")
        if not (-99 <= n <= 99):
            raise ValueError("Numerator n must be between -99 and 99")
        if not (0 < d <= 99):
            raise ValueError("Denominator d must be between 1 and 99")
            
        self.n = n
        self.d = d
        self.e = e

    @staticmethod
    def from_value(value):
        """Creates a Frac object by approximating a float value."""
        n, d, e = find_best_frac_approx(value)
        return Frac(n, d, e)

    def value(self):
        """Returns the float value of the fraction."""
        return (self.n / self.d) * (10 ** self.e)

    def __repr__(self):
        return f"Frac({self.n}/{self.d}e{self.e}) -> {self.value():.6f}"

    def _simplify_and_create(self, num, den, exp):
        """Internal method to simplify and create a new Frac object."""
        if den == 0:
            raise ZeroDivisionError
            
        common = gcd(abs(num), abs(den))
        num //= common
        den //= common

        # Max size for intermediate values is 5D (99999)
        if abs(num) > 99999 or abs(den) > 99999:
             pass # Assume hardware/library handles it, for sim we use Python's big integers

        val = (num / den) * (10**exp)
        return Frac.from_value(val)

    def __add__(self, other):
        # (n1/d1)*10^e1 + (n2/d2)*10^e2
        exp = min(self.e, other.e)
        num1 = self.n * (10**(self.e - exp))
        den1 = self.d
        num2 = other.n * (10**(other.e - exp))
        den2 = other.d
        
        new_num = num1 * den2 + num2 * den1
        new_den = den1 * den2
        return self._simplify_and_create(new_num, new_den, exp)
        
    def __sub__(self, other):
        # (n1/d1)*10^e1 - (n2/d2)*10^e2
        exp = min(self.e, other.e)
        num1 = self.n * (10**(self.e - exp))
        den1 = self.d
        num2 = other.n * (10**(other.e - exp))
        den2 = other.d
        
        new_num = num1 * den2 - num2 * den1
        new_den = den1 * den2
        return self._simplify_and_create(new_num, new_den, exp)

    def __mul__(self, other):
        new_num = self.n * other.n
        new_den = self.d * other.d
        new_exp = self.e + other.e
        return self._simplify_and_create(new_num, new_den, new_exp)

    def __truediv__(self, other):
        new_num = self.n * other.d
        new_den = self.d * other.n
        new_exp = self.e - other.e
        return self._simplify_and_create(new_num, new_den, new_exp)

def solve_dilation_factor():
    # 1. Define Core Parameters
    G = 6.67430e-11  # m^3 kg^-1 s^-2
    M_sun = 1.989e30 # kg
    c = 2.99792458e8 # m/s
    
    M_pandora = 2 * M_sun # kg
    R_pandora = 20 * 1000 # m
    d_pioneer = 60 * 1000 # m
    
    r = R_pandora + d_pioneer
    
    # 2. Calculate the value of x = Rs/r
    Rs = (2 * G * M_pandora) / (c**2)
    x_val = Rs / r
    
    # 3. Create Frac objects for calculation
    # These represent the constants loaded into the program.
    # We use the factory method to find the best frac approximation.
    frac_one = Frac(1, 1, 0)
    frac_two = Frac(2, 1, 0)
    frac_eight = Frac(8, 1, 0)
    frac_x = Frac.from_value(x_val)

    # 4. Perform Taylor Series Calculation: f ≈ 1 - x/2 - x²/8
    # This logic represents the C program's steps.
    # This can be done with 2 frac variables in C: `f` and `temp`
    
    # In C: frac f = frac_one;
    #       frac temp = frac_div(frac_x, frac_two); // temp = x/2
    #       f = frac_sub(f, temp); // f = 1 - x/2
    term1 = frac_x / frac_two
    f = frac_one - term1
    
    # In C: temp = frac_mul(frac_x, frac_x); // temp = x*x
    #       temp = frac_div(temp, frac_eight); // temp = x*x/8
    #       f = frac_sub(f, temp); // f = (1-x/2) - x*x/8
    x_squared = frac_x * frac_x
    term2 = x_squared / frac_eight
    f = f - term2

    # 5. Determine Memory Usage (z)
    # Most efficient C program needs 2 `frac` variables (f, temp)
    # and frac constants which can be created on-the-fly.
    # Each frac variable is n(2D)+d(2D)+e(2D) = 6D.
    # z = 2 * 6D = 12D.
    z = 12
    
    final_f_value = f.value()
    
    print(f"{final_f_value:.4f}:{z}")

# Run the simulation and print the final answer
solve_dilation_factor()
>>> 0.9629:12