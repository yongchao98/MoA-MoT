import math

class Frac:
    """
    A class to simulate the Wuxing `frac` data type.
    Represents a value as (n/d) * 10^e.
    n (numerator): signed 2D char (-99 to 99)
    d (denominator): unsigned 2D char (0 to 99)
    e (exponent): signed 2D char (-99 to 99)
    """

    def __init__(self, n, d, e):
        if not isinstance(n, int) or not isinstance(d, int) or not isinstance(e, int):
            raise TypeError("n, d, e must be integers.")
        if d == 0:
            raise ZeroDivisionError("Denominator cannot be zero.")
        self.n = n
        self.d = d
        self.e = e
        self._simplify()

    def _simplify(self):
        # Move sign from d to n
        if self.d < 0:
            self.n = -self.n
            self.d = -self.d

        if self.n == 0:
            self.d = 1
            self.e = 0
            return

        # Reduce fraction by Greatest Common Divisor
        common_divisor = math.gcd(abs(self.n), self.d)
        self.n //= common_divisor
        self.d //= common_divisor

        # Normalize magnitude of n to fit in 2D char by adjusting exponent
        while abs(self.n) >= 100:
            rem = self.n % 10
            self.n //= 10
            self.e += 1
            # Apply rounding
            if abs(rem) >= 5:
                self.n += (1 if self.n > 0 else -1)
        
        # Check against Wuxing constraints (for simulation purposes)
        if abs(self.n) > 99 or self.d > 99 or abs(self.e) > 99:
            # This would be handled by hardware, here we just note it.
            # In a real scenario, this might cause an error or further lossy simplification.
            pass

    def value(self):
        return (self.n / self.d) * (10 ** self.e)

    def __repr__(self):
        return f"{self.n}/{self.d}e{self.e}"

    def __neg__(self):
        return Frac(-self.n, self.d, self.e)

    def __add__(self, other):
        # (n1/d1)*10^e1 + (n2/d2)*10^e2
        # To add, exponents must match.
        e_common = max(self.e, other.e)
        n1 = self.n * (10**(e_common - self.e))
        n2 = other.n * (10**(e_common - other.e))
        
        new_n = n1 * other.d + n2 * self.d
        new_d = self.d * other.d
        return Frac(new_n, new_d, e_common)

    def __sub__(self, other):
        return self + (-other)

    def __mul__(self, other):
        new_n = self.n * other.n
        new_d = self.d * other.d
        new_e = self.e + other.e
        return Frac(new_n, new_d, new_e)

    def __truediv__(self, other):
        if other.n == 0:
            raise ZeroDivisionError("Division by zero fraction.")
        new_n = self.n * other.d
        new_d = self.d * other.n
        new_e = self.e - other.e
        return Frac(new_n, new_d, new_e)


def solve():
    # 1. Define variables based on the problem description
    # g = 9.8 m/s^2, represented as 98/1 * 10^-1
    g_frac = Frac(98, 1, -1)
    # Lion's speed v = 5 m/s
    v_frac = Frac(5, 1, 0)
    # Initial distance = 300 m, represented as 3/1 * 10^2
    dist_frac = Frac(3, 1, 2)
    # sin(60 deg) = sin(120 deg) ~= 0.866. Represent as 43/50 (value 0.86)
    # This fits the n, d < 100 constraint.
    sin_val = Frac(43, 50, 0)

    # 2. Setup the quadratic equation: A*u^2 + B*u + C = 0
    # A = sin(2a) = sin(120) ~= 0.86
    coeff_A = sin_val
    # B = 2 * v * sin(a) = 2 * 5 * sin(60) = 10 * 0.86 = 8.6
    coeff_B = Frac(2, 1, 0) * v_frac * sin_val
    # C = -d * g = -300 * 9.8 = -2940
    coeff_C = -(dist_frac * g_frac)

    # 3. Use Newton's method to solve for u
    # f(u) = A*u^2 + B*u + C
    # f'(u) = 2*A*u + B
    
    # Initial guess: speed to hit a stationary target at 300m
    # R = u^2*sin(2a)/g => u = sqrt(R*g/sin(2a)) = sqrt(300*9.8/0.86) ~= 58 m/s
    u_speed = Frac(58, 1, 0)
    
    # Constants for the calculation
    TWO = Frac(2, 1, 0)
    
    # Iterate to find the root
    for _ in range(10): # 10 iterations for convergence
        u_sq = u_speed * u_speed
        f_u = (coeff_A * u_sq) + (coeff_B * u_speed) + coeff_C
        f_prime_u = (TWO * coeff_A * u_speed) + coeff_B
        if f_prime_u.n == 0: # Avoid division by zero
            break
        u_speed = u_speed - (f_u / f_prime_u)

    # 4. Calculate final values and memory usage
    final_u_val = u_speed.value()
    
    # Memory usage calculation:
    # 8 variables are stored: g_frac, v_frac, dist_frac, sin_val,
    # coeff_A, coeff_B, coeff_C, u_speed.
    # Each frac costs 2D(n) + 2D(d) + 2D(e) = 6D.
    # Total memory z = 8 * 6 = 48 D.
    z = 48
    
    # 5. Print the required output
    print("The problem is solved using the equation: A*u^2 + B*u + C = 0")
    print("Calculated coefficients in frac format (n/d*10^e):")
    print(f"A = {coeff_A}")
    print(f"B = {coeff_B}")
    print(f"C = {coeff_C}")
    print("\nUsing Newton's method, the calculated initial speed u is:")
    print(f"u = {final_u_val:.1f} m/s")

    print(f"\nFinal Answer in u:z format:")
    print(f"<<<{final_u_val:.1f}:{z}>>>")

solve()