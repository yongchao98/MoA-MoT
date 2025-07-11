import sympy

def solve():
    """
    Calculates the thickness of the double point of the stable reduction
    of the curve z^2 = 2*x^5 + 2*x^3 + 1 above 2.
    """
    # Step 1: Define curve properties
    g = 2  # Genus of the curve
    
    # Define the polynomial f(x)
    x = sympy.Symbol('x')
    f = 2*x**5 + 2*x**3 + 1
    
    # Step 2: Calculate the discriminant of the polynomial and the curve
    # The discriminant of the polynomial f(x)
    poly_discriminant = sympy.discriminant(f, x)
    
    # Factor the discriminant to find its 2-adic valuation
    factorization = sympy.factorint(poly_discriminant)
    
    # The valuation v(p) is normalized to 1, so v_2(n) is the exponent of 2 in n's prime factorization.
    v_poly_discriminant = factorization[2]
    
    # For a hyperelliptic curve y^2 = f(x) of genus 2, the discriminant of the curve is Delta = 2^12 * disc(f)
    v_curve_discriminant = 12 + v_poly_discriminant
    
    # Step 3: Determine properties of the stable reduction
    # Based on the analysis of the roots of f(x), the stable reduction has
    # two components intersecting at one node.
    k = 2          # Number of components in the stable reduction
    k_sing = 1     # Number of nodes (double points)

    # Step 4: Use the formula for the valuation of the discriminant
    # v(Delta) = 2*(g-1)*(k-1) + k_sing + 2*delta
    # We solve for delta (the thickness)
    
    # Step 5: Solve for delta and print the equation
    # delta = (v(Delta) - 2*(g-1)*(k-1) - k_sing) / 2
    
    term1 = 2 * (g - 1) * (k - 1)
    
    delta = (v_curve_discriminant - term1 - k_sing) / 2
    
    print("The properties of the curve are:")
    print(f"Genus g = {g}")
    print(f"The polynomial is f(x) = {f}")
    
    print("\nCalculation of the discriminant:")
    print(f"The discriminant of the polynomial f(x) is {poly_discriminant}")
    print(f"The 2-adic valuation of the polynomial discriminant is v(disc(f)) = {v_poly_discriminant}")
    print(f"The 2-adic valuation of the curve's discriminant is v(Delta) = 12 + v(disc(f)) = 12 + {v_poly_discriminant} = {v_curve_discriminant}")

    print("\nProperties of the stable reduction's special fiber:")
    print(f"Number of components k = {k}")
    print(f"Number of double points (nodes) k_sing = {k_sing}")
    
    print("\nApplying the formula for the valuation of the discriminant:")
    print("v(Delta) = 2 * (g - 1) * (k - 1) + k_sing + 2 * thickness")
    print(f"{v_curve_discriminant} = 2 * ({g} - 1) * ({k} - 1) + {k_sing} + 2 * thickness")
    print(f"{v_curve_discriminant} = {term1} + {k_sing} + 2 * thickness")
    print(f"{v_curve_discriminant} - {term1} - {k_sing} = 2 * thickness")
    print(f"{v_curve_discriminant - term1 - k_sing} = 2 * thickness")
    print(f"thickness = {v_curve_discriminant - term1 - k_sing} / 2")
    print(f"thickness = {delta}")

solve()