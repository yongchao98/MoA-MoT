import math

class Frac:
    """
    Simulates the Wuxing 'frac' type: (n/d) * 10^e
    - n: signed char (2D, simulated with int)
    - d: unsigned char (2D, simulated with int)
    - e: signed char (2D, simulated with int)
    """
    def __init__(self, n, d, e):
        self.n = int(n)
        self.d = int(d)
        self.e = int(e)
        self._simplify()

    def _simplify(self):
        if self.d == 0:
            raise ZeroDivisionError("Denominator cannot be zero.")
        if self.n == 0:
            self.d = 1
            self.e = 0
            return
        
        common = math.gcd(self.n, self.d)
        self.n //= common
        self.d //= common

        if self.d < 0:
            self.n = -self.n
            self.d = -self.d
    
    def to_float(self):
        """Converts the fraction to a standard float for display."""
        return (self.n / self.d) * (10 ** self.e)

    def __add__(self, other):
        # Align exponents
        if self.e > other.e:
            e_new = self.e
            n2_adj = other.n * (10 ** (self.e - other.e))
            n1_adj = self.n
        else:
            e_new = other.e
            n1_adj = self.n * (10 ** (other.e - self.e))
            n2_adj = other.n

        n_new = n1_adj * other.d + n2_adj * self.d
        d_new = self.d * other.d
        return Frac(n_new, d_new, e_new)

    def __sub__(self, other):
        # Align exponents
        if self.e > other.e:
            e_new = self.e
            n2_adj = other.n * (10 ** (self.e - other.e))
            n1_adj = self.n
        else:
            e_new = other.e
            n1_adj = self.n * (10 ** (other.e - self.e))
            n2_adj = other.n

        n_new = n1_adj * other.d - n2_adj * self.d
        d_new = self.d * other.d
        return Frac(n_new, d_new, e_new)

    def __mul__(self, other):
        return Frac(self.n * other.n, self.d * other.d, self.e + other.e)

    def __truediv__(self, other):
        return Frac(self.n * other.d, self.d * other.n, self.e - other.e)

def solve_projectile_problem():
    # --- Variable Declaration and Memory Calculation ---
    # These variables are considered non-temporary for memory calculation.
    # 10 'frac' variables and 1 'int' variable.
    # sizeof(frac) = 6D, sizeof(int) = 5D
    # z = 10 * 6D + 1 * 5D = 65D
    
    # Problem constants
    mass = Frac(5, 10, 0)      # 0.5 kg (unused in kinematics but declared)
    g = Frac(98, 10, 0)       # 9.8 m/s^2
    v = Frac(5, 1, 0)         # 5 m/s
    dist = Frac(300, 1, 0)    # 300 m

    # Angle a = 60 degrees
    # cos(60) = 0.5
    cos60 = Frac(1, 2, 0)
    # sin(60) = sqrt(3)/2 ~= 0.866025. Use a precise fraction.
    sin60 = Frac(866, 1000, 0)

    # Helper constants
    two = Frac(2, 1, 0)

    # Binary search variables
    u_low = Frac(0, 1, 0)
    u_high = Frac(100, 1, 0)  # Initial guess for upper bound of speed u
    u_mid = Frac(0, 1, 0)

    # Loop counter
    iterations = 30 # An int variable
    for i in range(iterations):
        u_mid = (u_low + u_high) / two

        # Equation to satisfy: range_rock(u) - range_lion(u) = 0
        
        # t = (2 * u * sin(a)) / g
        t = (two * u_mid * sin60) / g
        
        # range_rock = u * cos(a) * t
        range_rock = u_mid * cos60 * t
        
        # range_lion = initial_dist - v * t
        range_lion = dist - (v * t)
        
        # error = range_rock - range_lion
        error = range_rock - range_lion
        
        if error.n > 0:
            u_high = u_mid
        else:
            u_low = u_mid
    
    # The final velocity is u_low (or u_high, they are very close)
    final_u = u_low

    # --- Outputting the results ---
    
    # Calculate final numbers for the equation as requested
    final_t = (two * final_u * sin60) / g
    final_range_rock = final_u * cos60 * final_t
    final_range_lion = dist - (v * final_t)
    
    print("Solved equation: rock_range = lion_position")
    print(f"The numbers in the final equation are:")
    print(f"rock_range = {final_range_rock.to_float():.4f}")
    print(f"lion_position = {final_range_lion.to_float():.4f}")
    print("-" * 20)

    # Calculate final answer u:z
    u_result = round(final_u.to_float(), 1)
    z_memory = 10 * 6 + 5 # 10 fracs, 1 int

    print(f"Final Answer (u:z):")
    print(f"{u_result}:{z_memory}")

# Execute the solution
solve_projectile_problem()