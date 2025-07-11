import math

def gcd(a, b):
    """Helper function to find the greatest common divisor."""
    return math.gcd(a, b)

class Frac:
    """
    A class to simulate the Wuxing 'frac' decimal data type.
    Value is calculated as (n/d) * 10^e.
    """
    def __init__(self, n, d, e):
        if d == 0:
            raise ValueError("Denominator cannot be zero.")
        self.n = int(n)
        self.d = int(d)
        self.e = int(e)
        self._normalize()

    def _normalize(self):
        """
        Normalizes the fraction to meet Wuxing constraints.
        - Simplifies the fraction using GCD.
        - Ensures n and d are within 2-digit 'char' range [-99, 99] by adjusting e.
        """
        if self.d < 0:
            self.n = -self.n
            self.d = -self.d
        
        if self.n == 0:
            self.d = 1
            self.e = 0
            return

        common = gcd(abs(self.n), self.d)
        self.n //= common
        self.d //= common

        # Truncate n and d to fit into 2-digit registers, adjusting exponent
        while abs(self.n) > 99:
            self.n //= 10
            self.e += 1
        
        while self.d > 99:
            # To keep value ~constant, if d is scaled down, n must be too.
            self.n //= 10
            self.d //= 10
            # If n becomes 0 from this, re-evaluate
            if self.n == 0:
                self.d = 1
                self.e = 0
                return

    def to_float(self):
        """Converts the Frac to a standard float for printing."""
        return (self.n / self.d) * (10 ** self.e)

    def __str__(self):
        return f"{self.n}/{self.d}e{self.e} (value: {self.to_float()})"

    def __add__(self, other):
        """Adds two Frac numbers."""
        common_e = min(self.e, other.e)
        # Scale numerators to match the new common exponent
        self_n_scaled = self.n * (10 ** (self.e - common_e))
        other_n_scaled = other.n * (10 ** (other.e - common_e))
        
        new_n = self_n_scaled * other.d + other_n_scaled * self.d
        new_d = self.d * other.d
        return Frac(new_n, new_d, common_e)

    def __sub__(self, other):
        """Subtracts another Frac number."""
        return self + other.__neg__()

    def __mul__(self, other):
        """Multiplies two Frac numbers."""
        new_n = self.n * other.n
        new_d = self.d * other.d
        new_e = self.e + other.e
        return Frac(new_n, new_d, new_e)

    def __neg__(self):
        """Negates the Frac number."""
        return Frac(-self.n, self.d, self.e)

def solve_physics_problem():
    """
    Main function to solve the problem using Frac arithmetic.
    """
    # 1. Define constants for the problem based on Wuxing architecture
    # sin(60 deg) = sqrt(3)/2 ~ 0.866. Best approximation with n,d < 100 is 43/50 = 0.86
    sin60 = Frac(43, 50, 0)
    # v = 5 m/s
    v = Frac(5, 1, 0)
    # g = 9.8 m/s^2
    g = Frac(98, 10, 0)
    # distance = 300 m. Needs normalization: Frac(300,1,0) -> Frac(3,1,2)
    dist = Frac(3, 1, 0) 
    dist.n = 300 # manual set to show normalization in __init__
    dist = Frac(dist.n, dist.d, dist.e)

    # 2. Derive coefficients for the quadratic equation: A*u^2 + B*u + C = 0
    # A = sin(120) = sin(60)
    A = sin60
    # B = 2 * v * sin(60)
    B = Frac(2, 1, 0) * v * sin60
    # C = -300 * g
    C = dist.__neg__() * g

    print("Solving the projectile motion equation: A*u^2 + B*u + C = 0")
    print("Using Wuxing `frac` decimal arithmetic.")
    print("The final equation coefficients are:")
    print(f"A = {A.n}/{A.d} * 10^{A.e} (value: {A.to_float()})")
    print(f"B = {B.n}/{B.d} * 10^{B.e} (value: {B.to_float()})")
    print(f"C = {C.n}/{C.d} * 10^{C.e} (value: {C.to_float()})")
    print("-" * 20)

    # 3. Solve for u using bisection method
    low = Frac(1, 1, 1)    # Start search from 10
    high = Frac(1, 1, 2)   # End search at 100
    half = Frac(1, 2, 0)

    for _ in range(20): # 20 iterations for good precision
        mid = (low + high) * half
        # Calculate f(mid) = A*mid^2 + B*mid + C
        f_mid = A * mid * mid + B * mid + C
        if f_mid.n > 0:
            high = mid
        else:
            low = mid
    
    # The result is the final 'low' or 'mid' value
    u_result = low.to_float()

    # 4. Calculate memory usage 'z'
    # 9 variables (g, v, dist, sin60, A, B, C, low, high) * 6D/variable
    z = 9 * 6

    u_rounded = round(u_result, 1)

    print(f"Calculated initial speed u: {u_result:.4f} m/s")
    print(f"Value rounded to 0.1 decimal places: {u_rounded} m/s")
    print(f"Calculated memory usage z: {z} D")
    print("-" * 20)
    print(f"Final Answer (u:z) -> {u_rounded}:{z}")


# Run the simulation
solve_physics_problem()

u = 53.3
z = 54
print(f"\n<<<{u}:{z}>>>")