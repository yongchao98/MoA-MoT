import sympy

def solve():
    """
    This program calculates the thickness of the double point of the stable reduction of the curve z^2 = 2*x^5 + 2*x^3 + 1.
    
    Plan:
    1. Define the polynomial f(x) from the curve equation.
    2. Calculate the discriminant of the polynomial f(x).
    3. Calculate the 2-adic valuation of the discriminant of f(x).
    4. Calculate the 2-adic valuation of the discriminant of the curve using the formula for hyperelliptic curves.
    5. Use the conductor formula v(Delta_C) = delta + epsilon to find the total thickness epsilon.
    6. Argue that delta = 0 because the special fiber is irreducible.
    7. Conclude that the thickness of the single double point is epsilon.
    """
    
    # Step 1: Define the polynomial
    x = sympy.Symbol('x')
    f = 2 * x**5 + 2 * x**3 + 1
    
    # Step 2: Calculate the discriminant of the polynomial
    # For a polynomial, the discriminant can be computed using sympy.discriminant
    disc_f = sympy.discriminant(f, x)
    
    # Step 3: Calculate the 2-adic valuation of the discriminant
    # The 2-adic valuation v(n) is the exponent of the highest power of 2 that divides n.
    n = int(disc_f)
    v_2_disc_f = 0
    if n == 0:
        # Valuation of 0 is infinite, but discriminant is non-zero here.
        v_2_disc_f = float('inf')
    else:
        while n % 2 == 0:
            v_2_disc_f += 1
            n //= 2
            
    # Step 4: Calculate the 2-adic valuation of the curve's discriminant
    # For a hyperelliptic curve y^2=f(x) of genus g, the discriminant is Delta_C = 2^(4(g-1)) * disc(f)
    # The degree of f is 5, so the genus g = floor((5-1)/2) = 2.
    g = 2
    # The valuation is v(Delta_C) = v(2^(4(g-1))) + v(disc(f))
    # With v(2)=1, v(2^k) = k.
    v_2_Delta_C = 4 * (g - 1) + v_2_disc_f
    
    # Step 5: Use the conductor formula v(Delta_C) = delta + epsilon
    # The analysis of the Newton polygon of f(x) shows all roots form a single cluster.
    # This implies the stable reduction has an irreducible special fiber.
    # For an irreducible special fiber, the number of components is 1.
    num_components = 1
    
    # Step 6: delta is the number of components - 1
    delta = num_components - 1
    
    # The total thickness epsilon is given by v(Delta_C) - delta
    epsilon = v_2_Delta_C - delta
    
    print(f"The polynomial is f(x) = {f}")
    print(f"The discriminant of f(x) is {disc_f}")
    print(f"The 2-adic valuation of the discriminant of f(x) is v(disc(f)) = {v_2_disc_f}")
    print(f"The genus of the curve is g = {g}")
    print(f"The 2-adic valuation of the discriminant of the curve is v(Delta_C) = 4*(g-1) + v(disc(f)) = 4*({g}-1) + {v_2_disc_f} = {v_2_Delta_C}")
    print(f"The special fiber of the stable reduction is irreducible, so the number of components is {num_components}.")
    print(f"This means delta = (number of components - 1) = {delta}")
    print(f"The conductor formula is v(Delta_C) = delta + epsilon, so {v_2_Delta_C} = {delta} + epsilon")
    print(f"The total thickness of singularities is epsilon = {epsilon}")
    print(f"Assuming there is a single double point, its thickness is {epsilon}.")
    
solve()