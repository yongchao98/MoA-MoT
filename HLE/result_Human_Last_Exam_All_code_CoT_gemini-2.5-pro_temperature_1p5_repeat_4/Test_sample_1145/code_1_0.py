import math

# --- Setup ---
# Quantum number
n = 2
# Interval limits as a fraction of the box length 'a'
x1_a = 0.495
x2_a = 0.505

# --- Calculation ---
# The probability P is the definite integral of |ψ(x)|² from x1 to x2.
# P = ∫ [from x1 to x2] (2/a) * sin²(n*π*x/a) dx
# The antiderivative of this function is:
# F(x) = (x/a) - (1/(2*n*π)) * sin(2*n*π*x/a)
# The probability is F(x2) - F(x1).
# Note that 'a' cancels out if we work with x_a = x/a.

# Let's define the antiderivative function in terms of z = x/a
def antiderivative(z, n_val):
    """Calculates the value of the integrated probability function."""
    return z - (1 / (2 * n_val * math.pi)) * math.sin(2 * n_val * math.pi * z)

# Evaluate the antiderivative at the interval boundaries
val_at_x2 = antiderivative(x2_a, n)
val_at_x1 = antiderivative(x1_a, n)

# The final probability is the difference
probability = val_at_x2 - val_at_x1

# --- Output ---
print("To find the probability that the particle is in the interval [0.495a, 0.505a] for n=2, we solve the following integral:")
print("P = ∫[from 0.495a to 0.505a] (2/a) * sin²(2*π*x/a) dx")
print("\nAfter integration, the expression to evaluate is:")
print("P = [x/a - sin(4*π*x/a) / (4*π)] evaluated from x=0.495a to x=0.505a")
print("\nThis results in the final equation:")
# We reconstruct the equation with all the numbers for clarity, as requested.
# Equation form: P = (x2/a - x1/a) - (1/(4*π)) * [sin(4*π*x2/a) - sin(4*π*x1/a)]
sin_term_2_val = math.sin(2 * n * math.pi * x2_a)
sin_term_1_val = math.sin(2 * n * math.pi * x1_a)
coeff_val = 1 / (2 * n * math.pi)

print(f"P = ({x2_a} - {x1_a}) - ({round(coeff_val, 5)}) * [sin(4*π*{x2_a}) - sin(4*π*{x1_a})]")
print("\nCalculating the numerical values of each part:")
print(f"sin(4*π*{x2_a}) = sin({round(4*x2_a, 3)}π) = {sin_term_2_val}")
print(f"sin(4*π*{x1_a}) = sin({round(4*x1_a, 3)}π) = {sin_term_1_val}")
print(f"Difference in sines = {sin_term_2_val} - ({sin_term_1_val}) = {sin_term_2_val - sin_term_1_val}")

print("\nPlugging these into the equation:")
print(f"P = {x2_a - x1_a} - {coeff_val} * ({sin_term_2_val - sin_term_1_val})")
print(f"P = {x2_a - x1_a} - {coeff_val * (sin_term_2_val - sin_term_1_val)}")

print(f"\nThe final probability is:")
print(probability)

<<<6.579737233890388e-06>>>