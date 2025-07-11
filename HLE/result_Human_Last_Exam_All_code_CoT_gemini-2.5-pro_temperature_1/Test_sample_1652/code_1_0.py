import math

def gcd(a, b):
    """Computes the greatest common divisor of a and b."""
    return math.gcd(a, b)

class Frac:
    """
    Simulates the Wuxing 'frac' data type.
    Value is (n/d) * 10^e.
    n: signed char (-127 to 127)
    d: unsigned char (1 to 255)
    e: signed char (-127 to 127)
    """
    def __init__(self, n, d, e):
        if d == 0:
            raise ZeroDivisionError("Denominator cannot be zero.")
        if d < 0: # Keep denominator positive
            n = -n
            d = -d
        self.n = n
        self.d = d
        self.e = e
        self.normalize()

    def normalize(self):
        """Simplifies the fraction and adjusts n, d, e to fit char limits."""
        if self.n == 0:
            self.d = 1
            self.e = 0
            return

        # Simplify n/d using GCD
        common = gcd(abs(self.n), self.d)
        self.n //= common
        self.d //= common

        # Adjust exponent to bring n and d into their range
        while abs(self.n) > 127 or self.d > 255:
            self.n = round(self.n / 10)
            self.d = round(self.d / 10)
            if self.d == 0: self.d = 1
            self.e += 1
        
        while self.e > 0 and abs(self.n * 10) <= 127 and self.d * 10 <= 255:
            # this check is tricky, for now we don't expand
            # as it might lose precision from the other side.
            # Let's prioritize fitting into the char types.
            break

    def value(self):
        """Returns the float value of the fraction."""
        return (self.n / self.d) * (10 ** self.e)

    def __repr__(self):
        return f"{self.n}/{self.d}e{self.e}"

    def __add__(self, other):
        # To add (n1/d1)*10^e1 + (n2/d2)*10^e2, align exponents
        e_diff = self.e - other.e
        if e_diff > 0:
            n1 = self.n * (10**e_diff)
            d1 = self.d
            n2, d2 = other.n, other.d
        else:
            n1 = self.n
            d1 = self.d
            n2 = other.n * (10**(-e_diff))
            d2 = other.d
        
        new_n = n1 * d2 + n2 * d1
        new_d = d1 * d2
        return Frac(new_n, new_d, min(self.e, other.e))

    def __sub__(self, other):
        return self + Frac(-other.n, other.d, other.e)

    def __mul__(self, other):
        return Frac(self.n * other.n, self.d * other.d, self.e + other.e)

    def __truediv__(self, other):
        return Frac(self.n * other.d, self.d * other.n, self.e - other.e)

def frac_sqrt(f, iterations=20):
    """Calculates square root of a Frac using the Babylonian method."""
    if f.n < 0:
        raise ValueError("Cannot calculate square root of a negative number.")
    # Initial guess
    val = f.value()
    guess_val = math.sqrt(val)
    # Create a Frac from the guess value. Crude but effective for simulation.
    e = math.floor(math.log10(guess_val)) if guess_val > 0 else 0
    n = round(guess_val / (10**e) * 100)
    d = 100
    g = Frac(n, d, e)
    
    two = Frac(2, 1, 0)
    for _ in range(iterations):
        g = (g + f / g) / two
    return g

# --- Main Program ---

# Define constants for the Wuxing system
v = Frac(5, 1, 0)       # 5 m/s
D = Frac(3, 1, 2)       # 300 m
g = Frac(98, 10, 0)     # 9.8 m/s^2

# Use the derived formula: u = (-v + sqrt(v^2 + 2Dg)) / sqrt(3)
# Calculate terms for the equation
two = Frac(2, 1, 0)
three = Frac(3, 1, 0)
neg_v = Frac(-v.n, v.d, v.e)

term_in_sqrt = v*v + two*D*g
sqrt_of_term = frac_sqrt(term_in_sqrt)
sqrt_of_3 = frac_sqrt(three)

numerator = neg_v + sqrt_of_term
u = numerator / sqrt_of_3

# Calculate memory usage
z = 24 # 4 Frac variables (v, D, g, u) * 6D/Frac = 24D

# Print the final equation with calculated values
print(f"Final Equation: u = ( {neg_v.value()} + sqrt({term_in_sqrt.value()}) ) / {sqrt_of_3.value()}")
print(f"Calculated u = ( {neg_v.value()} + {sqrt_of_term.value()} ) / {sqrt_of_3.value()}")
print(f"Result: u = {u.value():.1f} m/s, Memory Usage z = {z} D")

u_rounded = round(u.value(), 1)
final_answer = f"{u_rounded}:{z}"
print(f"<<<{final_answer}>>>")