import sympy

def get_valuation(n, p):
    """Calculates the p-adic valuation of an integer n."""
    if n == 0:
        return float('inf')
    val = 0
    n = abs(n)
    while n % p == 0:
        val += 1
        n //= p
    return val

# Define the variable and the polynomial
x = sympy.Symbol('x')
f = 2*x**5 + 2*x**3 + 1

# Calculate the discriminant
disc_f = int(sympy.discriminant(f, x))

# Calculate the 2-adic valuation of the discriminant
v_disc = get_valuation(disc_f, 2)

# Sizes of the root clusters from the Newton polygon
n1 = 3
n2 = 2

# Calculate the conductor exponent (thickness)
f = v_disc - ((n1 - 1) + (n2 - 1))

# Output the final equation
print(f"The thickness is calculated by the formula: f = v(Disc(f)) - (n1 - 1) - (n2 - 1)")
print(f"v(Disc(f)) = v({disc_f}) = {v_disc}")
print(f"n1 = {n1}")
print(f"n2 = {n2}")
print(f"f = {v_disc} - ({n1} - 1) - ({n2} - 1) = {v_disc} - {n1-1} - {n2-1} = {f}")