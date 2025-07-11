import math

class TitanNumber:
    """Represents a number in Titan's format: (n/d) * 10^e."""
    def __init__(self, n, d=1, e=0):
        # Initial check for direct creation
        if not (0 <= n <= 15 and 1 <= d <= 15):
             # This constructor is also used internally before reduction
             pass
        self.n = n
        self.d = d
        self.e = e

    def to_float(self):
        return (self.n / self.d) * (10 ** self.e)

    def __str__(self):
        return f"({self.n}/{self.d})e{self.e}"

# Precompute valid fractions for fast lookup
VALID_FRACTIONS = []
for d in range(1, 16):
    for n in range(0, 16):
        VALID_FRACTIONS.append(((n / d), n, d))
# Sort by value, then by denominator to break ties consistently
VALID_FRACTIONS.sort(key=lambda x: (x[0], x[2]))

def find_best_fraction(value):
    """Finds the best n/d representation for a value."""
    best_n, best_d = 0, 1
    min_err = float('inf')
    for val, n, d in VALID_FRACTIONS:
        err = abs(val - value)
        if err < min_err:
            min_err = err
            best_n, best_d = n, d
    return best_n, best_d

def reduce_value(value):
    """Reduces a float value to the best TitanNumber representation."""
    if value == 0:
        return TitanNumber(0, 1, 0)
    
    # Determine the best exponent by checking a small range around the log10 value
    # to find the mantissa that has the best fractional approximation.
    best_n, best_d, best_e = 0, 1, 0
    min_overall_err = float('inf')

    # The ideal exponent
    ideal_exp = math.floor(math.log10(abs(value))) if value != 0 else 0
    
    # Search exponents around the ideal one for a better fit
    for exp_offset in range(-1, 2):
      current_exp = ideal_exp + exp_offset
      mantissa = value / (10**current_exp)
      
      n, d = find_best_fraction(mantissa)
      
      err = abs(value - (n/d) * 10**current_exp)
      if err < min_overall_err:
        min_overall_err = err
        best_n, best_d, best_e = n, d, current_exp
        
    return TitanNumber(best_n, best_d, best_e)

def simplify_and_reduce(n, d, e):
    """Simplifies a raw fraction and reduces if it violates the 4-bit constraint."""
    if n == 0:
        return TitanNumber(0, 1, 0)
    common_divisor = math.gcd(n, d)
    n //= common_divisor
    d //= common_divisor

    if n > 15 or d > 15:
        return reduce_value((n / d) * (10 ** e))
    else:
        return TitanNumber(n, d, e)

def mul_titan(t1, t2):
    n = t1.n * t2.n
    d = t1.d * t2.d
    e = t1.e + t2.e
    return simplify_and_reduce(n, d, e)

def div_titan(t1, t2):
    if t2.n == 0:
        raise ZeroDivisionError
    n = t1.n * t2.d
    d = t1.d * t2.n
    e = t1.e - t2.e
    return simplify_and_reduce(n, d, e)

def add_titan_simplified(t1, t2):
    val1 = t1.to_float()
    val2 = t2.to_float()
    # Per rule, drop negligible terms. Define negligible as < 0.1% of the other.
    if abs(val1) < abs(val2) * 0.001:
        return t2
    if abs(val2) < abs(val1) * 0.001:
        return t1
    return reduce_value(val1 + val2)

def sqrt_titan(A, initial_guess):
    x = initial_guess
    one_half = TitanNumber(1, 2)
    for _ in range(4): # Iterate to find the solution
        A_div_x = div_titan(A, x)
        sum_val = add_titan_simplified(x, A_div_x)
        x_new = mul_titan(one_half, sum_val)
        if x_new.to_float() == x.to_float(): break
        x = x_new
    return x

# --- Main Calculation ---
G_true = 6.6743e-11
PI_true = 3.14159265
r_core_true, R_true = 1e5, 2e6
rho_core_true, rho_shell_true = 1200, 300
M_total_true = (4/3)*PI_true * (rho_core_true*r_core_true**3 + rho_shell_true*(R_true**3 - r_core_true**3))
v_e_true = (2 * G_true * M_total_true / R_true)**0.5

print("Starting Titan escape velocity calculation for Pandora.\n")
print("v_e^2 = (8 * pi * G / 3) * [ (rho_core - rho_shell) * r_core^3 / R + rho_shell * R^2 ]\n")

# Step 1: Define constants and parameters with best possible approximations
G = TitanNumber(13, 2, -11) # 6.5e-11
PI = TitanNumber(13, 4)     # 3.25
delta_rho = TitanNumber(9, 1, 2)  # 900
rho_shell = TitanNumber(3, 1, 2)  # 300
r_core = TitanNumber(1, 1, 5)     # 1e5
R = TitanNumber(2, 1, 6)          # 2e6

# Step 2: Calculate the two terms in the sum
r_core_cubed = mul_titan(mul_titan(r_core, r_core), r_core)
term_A_num = mul_titan(delta_rho, r_core_cubed)
Term_A = div_titan(term_A_num, R)

R_squared = mul_titan(R, R)
Term_B = mul_titan(rho_shell, R_squared)

# Step 3: Sum the terms (Term A is negligible and dropped)
Bracket_Sum = add_titan_simplified(Term_A, Term_B)

# Step 4: Calculate the constant C = 8*pi*G/3
eight_pi = mul_titan(TitanNumber(8, 1), PI)
eight_pi_G = mul_titan(eight_pi, G)
C = div_titan(eight_pi_G, TitanNumber(3, 1))

# Step 5: Calculate v_e^2
v_e_squared = mul_titan(C, Bracket_Sum)

# Step 6: Calculate v_e using Newton-Raphson method
initial_guess = TitanNumber(8, 1, 2) # Guess = 800 m/s
v_e_titan = sqrt_titan(v_e_squared, initial_guess)
v_e_titan_val = v_e_titan.to_float()

error = abs(v_e_titan_val - v_e_true)

# Final Output
print("--- Final Equation Numbers ---")
print(f"v_e^2 = C * (Sum of terms)")
print(f"C (8*pi*G/3) calculated as: {C.n}/{C.d}e{C.e}")
print(f"Sum of terms (dominated by shell) calculated as: {Bracket_Sum.n}/{Bracket_Sum.d}e{Bracket_Sum.e}")
print(f"\nFinal calculated v_e^2 = {v_e_squared.n}/{v_e_squared.d}e{v_e_squared.e}")
print(f"Resulting v_e from sqrt = {v_e_titan.n}/{v_e_titan.d}e{v_e_titan.e} which is {v_e_titan_val:.2f} m/s")
print(f"True value is {v_e_true:.2f} m/s.")
print(f"Final error is {error:.2f} m/s.")

print(f"\n\n<<<Y{error:.2f}>>>")