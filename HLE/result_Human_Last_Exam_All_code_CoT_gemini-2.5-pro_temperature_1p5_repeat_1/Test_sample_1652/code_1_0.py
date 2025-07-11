import math

class Frac:
    """
    A class to simulate the Wuxing 'frac' data type.
    It handles rational number arithmetic and simplification.
    The 'e' (exponent) is omitted as it can be folded into the numerator/denominator.
    """
    def __init__(self, n, d=1):
        if d == 0:
            raise ZeroDivisionError
        # Simplify the fraction using the greatest common divisor
        common = math.gcd(int(n), int(d))
        self.n = int(n // common)
        self.d = int(d // common)
        # Conventionally keep the denominator positive
        if self.d < 0:
            self.n = -self.n
            self.d = -self.d

    def __repr__(self):
        return f"{self.n}/{self.d}"

    def value(self):
        """Returns the floating-point value of the fraction."""
        return self.n / self.d

    def __add__(self, other):
        if not isinstance(other, Frac): other = Frac(other)
        new_n = self.n * other.d + other.n * self.d
        new_d = self.d * other.d
        return Frac(new_n, new_d)

    def __sub__(self, other):
        if not isinstance(other, Frac): other = Frac(other)
        new_n = self.n * other.d - other.n * self.d
        new_d = self.d * other.d
        return Frac(new_n, new_d)

    def __mul__(self, other):
        if not isinstance(other, Frac): other = Frac(other)
        new_n = self.n * other.n
        new_d = self.d * other.d
        return Frac(new_n, new_d)

    def __truediv__(self, other):
        if not isinstance(other, Frac): other = Frac(other)
        new_n = self.n * other.d
        new_d = self.d * other.n
        return Frac(new_n, new_d)

def solve_projectile_speed():
    """
    Calculates the initial speed 'u' and memory usage 'z'.
    """
    # 1. Define physics and problem constants as fractions
    g = Frac(98, 10)         # Gravitational acceleration, 9.8 m/s^2
    v = Frac(5)              # Lion's speed, 5 m/s
    distance = Frac(300)     # Initial distance, 300 m
    
    # Angle a = 60 degrees. We need its sin and cos.
    cos60 = Frac(1, 2)
    # sin(60) = sqrt(3)/2. We approximate sqrt(3) ~= 1.732.
    # So, sin(60) ~= 1.732 / 2 = 0.866 = 866/1000
    sin60 = Frac(866, 1000)

    # 2. Formulate the quadratic equation A*u^2 + B*u - C = 0
    # A = (2 * sin(a) * cos(a)) / g
    A = (Frac(2) * sin60 * cos60) / g
    # B = (2 * v * sin(a)) / g
    B = (Frac(2) * v * sin60) / g
    # C is the distance
    C = distance

    # 3. Solve for 'u' using the bisection method
    u_low = Frac(0)
    u_high = Frac(100) # Initial guess for the upper bound of the speed

    for _ in range(50): # 50 iterations provide sufficient precision
        u_mid = (u_low + u_high) / Frac(2)
        # Evaluate f(u) = A*u^2 + B*u - C
        f_mid_val = (A * u_mid * u_mid) + (B * u_mid) - C
        
        if f_mid_val.value() < 0:
            u_low = u_mid
        else:
            u_high = u_mid
            
    # The result is the midpoint of the final interval
    result_u_frac = (u_low + u_high) / Frac(2)
    result_u = result_u_frac.value()

    # 4. Calculate memory usage 'z'
    # The main non-temporary variables are: g, v, distance, sin60, cos60, A, B, u_low, u_high
    num_frac_variables = 9
    # Each 'frac' type costs 6D (2D for n, 2D for d, 2D for e)
    size_per_frac = 6
    z = num_frac_variables * size_per_frac

    # 5. Print the output as specified
    print(f"Equation: ({A}) * u^2 + ({B}) * u - ({C}) = 0")
    print("Result u:z")
    print(f"{result_u:.1f}:{z}")


solve_projectile_speed()
<<<53.5:54>>>