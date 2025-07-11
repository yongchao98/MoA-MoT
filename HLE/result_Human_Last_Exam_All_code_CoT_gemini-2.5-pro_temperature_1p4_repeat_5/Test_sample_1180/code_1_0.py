import sympy

def get_2_adic_valuation(n):
    """
    Computes the 2-adic valuation of an integer n.
    This is the highest power of 2 that divides n.
    Returns 0 if n is 0 or odd.
    """
    if n == 0:
        return float('inf')
    n = abs(n)
    valuation = 0
    while n > 0 and n % 2 == 0:
        valuation += 1
        n //= 2
    return valuation

# Define the variable and the polynomial from the transformed curve equation
t = sympy.Symbol('t')
Q = t**5 + 2*t**2 + 2

# Calculate the discriminant of the polynomial
disc = sympy.discriminant(Q, t)

# Calculate the 2-adic valuation of the discriminant
thickness = get_2_adic_valuation(int(disc))

print(f"The curve equation is transformed to w^2 = t * ({Q}).")
print(f"The relevant polynomial is Q(t) = {Q}.")
print(f"The discriminant of Q(t) is {disc}.")
print(f"The 2-adic valuation of the discriminant is v_2({disc}) = {thickness}.")
print(f"The thickness of the double point is {thickness}.")