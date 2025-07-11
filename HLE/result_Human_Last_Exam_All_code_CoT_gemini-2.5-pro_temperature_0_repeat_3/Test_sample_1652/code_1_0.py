import math

def gcd(a, b):
    """Computes the greatest common divisor of a and b."""
    return math.gcd(a, b)

class Frac:
    """
    A class to represent fractions and perform arithmetic, simulating
    the Wuxing architecture's 'frac' type. It avoids floating-point numbers.
    """
    def __init__(self, n, d=1):
        if d == 0:
            raise ZeroDivisionError
        common = gcd(n, d)
        self.n = n // common
        self.d = d // common

    def __repr__(self):
        return f"({self.n}/{self.d})"

    def to_float(self):
        return self.n / self.d

    def __add__(self, other):
        new_n = self.n * other.d + other.n * self.d
        new_d = self.d * other.d
        return Frac(new_n, new_d)

    def __sub__(self, other):
        new_n = self.n * other.d - other.n * self.d
        new_d = self.d * other.d
        return Frac(new_n, new_d)

    def __mul__(self, other):
        new_n = self.n * other.n
        new_d = self.d * other.d
        return Frac(new_n, new_d)

    def __truediv__(self, other):
        new_n = self.n * other.d
        new_d = self.d * other.n
        return Frac(new_n, new_d)

def solve_monkey_problem():
    """
    Solves the physics problem using the Frac class to adhere to
    Wuxing architecture constraints.
    """
    # Define constants as integers or fractions
    # g = 9.8 m/s^2
    g = Frac(98, 10)
    # v = 5 m/s
    v = Frac(5, 1)
    # d = 300 m
    d = Frac(300, 1)
    # a = 60 degrees. We need sin(60) and cos(60).
    # sin(60) = sqrt(3)/2 ~ 0.866
    sin_a = Frac(866, 1000)
    # cos(60) = 1/2
    cos_a = Frac(1, 2)

    # Memory usage calculation (z)
    # frac g (6D), int v (5D), int d (5D), frac sin_a (6D), frac cos_a (6D), frac u (6D)
    # Note: v and d are conceptually integers but are converted to Frac for calculation.
    # We count the variables needed in the C-style program.
    z = 6 + 5 + 5 + 6 + 6 + 6

    # Iterative search for the best initial speed u with 0.1 precision
    best_u = Frac(0, 1)
    min_diff = float('inf')

    # Search u from 50.0 to 60.0, as this is a reasonable range
    for u_scaled in range(500, 600):
        u = Frac(u_scaled, 10)

        # The physics equation is: u^2*sin(2a) + 2*u*v*sin(a) - d*g = 0
        # where sin(2a) = 2*sin(a)*cos(a)
        
        term1 = u * u * Frac(2,1) * sin_a * cos_a
        term2 = Frac(2,1) * u * v * sin_a
        term3 = d * g
        
        current_diff = term1 + term2 - term3
        
        if abs(current_diff.to_float()) < min_diff:
            min_diff = abs(current_diff.to_float())
            best_u = u

    # The final equation is A*u^2 + B*u - C = 0
    # A = 2*sin(a)*cos(a)
    # B = 2*v*sin(a)
    # C = d*g
    A = Frac(2,1) * sin_a * cos_a
    B = Frac(2,1) * v * sin_a
    C = d * g
    
    u_val = best_u.to_float()
    
    print("Final Equation with solved u:")
    print(f"{A} * {u_val:.1f}^2 + {B} * {u_val:.1f} - {C} = {min_diff:.4f} (closest to 0)")
    
    # Final answer in the format u:z
    print("\nFinal Answer:")
    print(f"u:{u_val:.1f}")
    print(f"z:{z}")
    
    final_answer_string = f"<<<{u_val:.1f}:{z}>>>"
    # This print is for the final wrapped answer as requested.
    print(final_answer_string)


solve_monkey_problem()