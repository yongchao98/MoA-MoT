import math

def gcd(a, b):
    """Computes the greatest common divisor of two integers."""
    return math.gcd(a, b)

class Frac:
    """
    A class to represent fractions and perform arithmetic under Titan's 5-bit constraints.
    Numerators and denominators must be integers between 0 and 31.
    """
    def __init__(self, n, d):
        if not (0 <= n <= 31 and 0 <= d <= 31):
            raise ValueError(f"Numerator {n} or Denominator {d} out of 5-bit range (0-31)")
        self.n = n
        self.d = d

    def __mul__(self, other):
        """
        Multiplies two fractions with cross-cancellation to prevent intermediate overflow.
        Checks if the final result still fits within the 5-bit constraints.
        """
        n1, d1 = self.n, self.d
        n2, d2 = other.n, other.d

        # Cross-simplify before multiplying to avoid large intermediate numbers
        c1 = gcd(n1, d2)
        n1 //= c1
        d2 //= c1

        c2 = gcd(n2, d1)
        n2 //= c2
        d1 //= c2
        
        # Perform the final multiplication
        final_n = n1 * n2
        final_d = d1 * d2

        return Frac(final_n, final_d)

    def __repr__(self):
        return f"{self.n}/{self.d}"

# --- Main Titan Calculation ---

print("Step 1: Calculate the mass of the rock.")

# Constants from the problem description, represented as Titan fractions
# m = (4/3) * pi * r^3 * rho
# m = (4/3) * pi * (1/8) * (9/10)
# First, calculate the rational part of the mass formula: (4/3) * (1/8) * (9/10)
term1 = Frac(4, 3)
term2 = Frac(1, 8)
term3 = Frac(9, 10)

m_rational_part = term1 * term2 * term3
print(f"Rational part of mass calculation: ({term1}) * ({term2}) * ({term3}) = {m_rational_part}")
# So, mass m = (3/20) * pi

# Choose a Titan-compatible approximation for pi
pi_approx = Frac(26, 9)
print(f"Using approximation pi ≈ {pi_approx}")

# Calculate mass m
m = m_rational_part * pi_approx
print(f"Calculated mass m = ({m_rational_part}) * ({pi_approx}) = {m}\n")

print("Step 2: Calculate the required force F.")
# The formula for y=10m is F = 2 * m * g * sqrt(2)

# Choose Titan-compatible approximations for g and sqrt(2)
g_approx = Frac(10, 1)
sqrt2_approx = Frac(3, 2)
print(f"Using approximation g ≈ {g_approx}")
print(f"Using approximation sqrt(2) ≈ {sqrt2_approx}")

# Group calculations carefully to avoid overflow
# F = (2 * g) * (m * sqrt(2))
term_A = Frac(2, 1) * g_approx
print(f"Intermediate term (2 * g) = (2/1) * ({g_approx}) = {term_A}")
term_B = m * sqrt2_approx
print(f"Intermediate term (m * sqrt(2)) = ({m}) * ({sqrt2_approx}) = {term_B}")

# Final force calculation
F_titan = term_A * term_B
print("\nFinal force equation on Titan:")
print(f"F = ({term_A.n}/{term_A.d}) * ({term_B.n}/{term_B.d}) = {F_titan.n}/{F_titan.d}")

# Step 3: Calculate the error of the Titan result.
F_titan_val = F_titan.n / F_titan.d

# Calculate the true value using high-precision constants
pi_true = 3.1415926535
m_true = 0.15 * pi_true
g_true = 9.8
sqrt2_true = 1.4142135623
F_true = 2 * m_true * g_true * sqrt2_true
error = abs(F_titan_val - F_true)

print(f"\nTitan calculated Force = {F_titan_val:.3f} N")
print(f"True physical Force    = {F_true:.3f} N")
print(f"Absolute error e       = {error:.3f}")

# Final Answer Determination
# The calculated force F=13N falls within the required range to hit the coconut [12.94N, 13.20N]
# The calculation is possible. The error 'e' is the absolute difference to the central target's true value.
final_answer = f"Y[{error:.3f}]"
print(f"\nFinal Answer: {final_answer}")