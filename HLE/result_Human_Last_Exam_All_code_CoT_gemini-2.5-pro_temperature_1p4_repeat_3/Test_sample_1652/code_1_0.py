import math

# Simulate the Wuxing 'frac' type using a Python class.
# The Wuxing architecture's built-in simplification and overflow protection
# are mimicked by using Python's arbitrary-precision integers and a gcd helper.
class Fraction:
    def __init__(self, n, d=1):
        if d == 0:
            raise ZeroDivisionError
        common = math.gcd(n, d)
        self.n = n // common
        self.d = d // common
    
    def __repr__(self):
        return f"{self.n}/{self.d}"

    def __add__(self, other):
        return Fraction(self.n * other.d + other.n * self.d, self.d * other.d)

    def __sub__(self, other):
        return Fraction(self.n * other.d - other.n * self.d, self.d * other.d)

    def __mul__(self, other):
        return Fraction(self.n * other.n, self.d * other.d)

    def __truediv__(self, other):
        return Fraction(self.n * other.d, self.d * other.n)

    def to_float(self):
        return self.n / self.d

def solve_physics_problem():
    """
    Solves for the initial speed u and calculates memory usage z.
    """
    # Define constants based on the problem description.
    # We use fraction approximations suitable for the Wuxing architecture.
    v = Fraction(5)         # Lion's speed in m/s
    distance = Fraction(300)  # Initial distance in m
    g = Fraction(98, 10)      # Gravity g = 9.8 m/s^2, simplified to 49/5
    # angle a = 60 degrees. sin(60) ~ 0.866, cos(60) = 0.5
    sin60 = Fraction(866, 1000) # Simplified to 433/500
    cos60 = Fraction(1, 2)
    two = Fraction(2)
    
    # The equation to solve for u is: t_flight = t_intercept
    # (2 * u * sin(60)) / g = distance / (u * cos(60) + v)
    print("Equation to solve for u: (2 * u * sin(60)) / g = distance / (u * cos(60) + v)")
    print("\nUsing the following Wuxing frac representations:")
    print(f"v = {v}")
    print(f"distance = {distance}")
    print(f"g = {g}")
    print(f"sin(60) = {sin60}")
    print(f"cos(60) = {cos60}")
    
    # Binary search for u, as sqrt() is not available to solve the quadratic form directly.
    low_u = Fraction(0)
    high_u = Fraction(1000) # A sufficiently large upper bound for speed

    for _ in range(100): # 100 iterations for high precision
        mid_u = (low_u + high_u) / two
        
        # t_flight = (2 * u * sin(60)) / g
        t_flight_num = two * mid_u * sin60
        t_flight = t_flight_num / g
        
        # t_intercept = distance / (u * cos(60) + v)
        t_intercept_den = mid_u * cos60 + v
        t_intercept = distance / t_intercept_den

        if t_flight.to_float() < t_intercept.to_float():
            low_u = mid_u
        else:
            high_u = mid_u

    u_result = high_u.to_float()

    # Calculate memory usage (z)
    # Variables required: v (int), distance (int), g (frac), sin60 (frac), cos60 (frac)
    # Sizes: int = 5D, frac = 2D(n) + 2D(d) + 2D(e) = 6D
    mem_v = 5       # stored as int
    mem_distance = 5 # stored as int
    mem_g = 6        # must be a frac
    mem_sin60 = 6    # must be a frac
    mem_cos60 = 6    # must be a frac
    z_usage = mem_v + mem_distance + mem_g + mem_sin60 + mem_cos60
    
    u_rounded = round(u_result, 1)

    print(f"\nCalculated initial speed u (rounded to 0.1): {u_rounded}")
    print(f"Memory usage for variables z (in D): {z_usage}")
    
    print(f"\nFinal answer: {u_rounded}:{z_usage}")
    

solve_physics_problem()
<<<53.5:28>>>