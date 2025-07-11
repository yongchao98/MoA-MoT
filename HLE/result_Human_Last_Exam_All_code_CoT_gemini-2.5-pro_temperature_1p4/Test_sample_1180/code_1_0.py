import sympy

def solve_curve_thickness():
    """
    Calculates the thickness of the double point of the stable reduction of the given curve.
    """
    x = sympy.Symbol('x')
    f = 2*x**5 + 2*x**3 + 1
    g = 2

    print(f"The curve is defined by z^2 = f(x), where f(x) = {f}.")
    print(f"This is a hyperelliptic curve of genus g = {g}.")
    print("-" * 20)

    # Step 1: Calculate the discriminant of the polynomial f(x)
    delta_f = sympy.discriminant(f, x)
    print(f"Step 1: The discriminant of the polynomial f(x) is Delta_f = {delta_f}.")

    # Helper function for 2-adic valuation
    def val_2(n):
        if n == 0:
            return float('inf')
        n = abs(int(n))
        count = 0
        while n > 0 and n % 2 == 0:
            count += 1
            n //= 2
        return count

    v2_delta_f = val_2(delta_f)
    print(f"The 2-adic valuation of Delta_f is v2(Delta_f) = {v2_delta_f}.")
    print("-" * 20)
    
    # Step 2: Calculate the 2-adic valuation of the curve's discriminant
    # Formula: v2(Delta_C) = 4*(g-1) + v2(Delta_f)
    v2_delta_c_val = 4 * (g - 1) + v2_delta_f
    
    print("Step 2: The 2-adic valuation of the curve's minimal discriminant, v2(Delta_C), is calculated.")
    print(f"Using the formula v2(Delta_C) = 4*(g-1) + v2(Delta_f):")
    print(f"v2(Delta_C) = 4*({g}-1) + {v2_delta_f} = {v2_delta_c_val}")
    print("-" * 20)

    # Step 3: Apply Saito's formula for curves with potentially good reduction
    # v2(Delta_C) = 2g - 2 + delta_tot
    print("Step 3: Apply Saito's formula for a curve of genus g with potentially good reduction:")
    print("v2(Delta_C) = 2*g - 2 + delta_tot, where delta_tot is the total thickness of the nodes.")
    
    # Step 4: Solve for delta_tot
    rhs_const = 2*g - 2
    delta_tot = v2_delta_c_val - rhs_const

    print("\nSolving for delta_tot:")
    print(f"{v2_delta_c_val} = 2*{g} - 2 + delta_tot")
    print(f"{v2_delta_c_val} = {rhs_const} + delta_tot")
    print(f"delta_tot = {v2_delta_c_val} - {rhs_const} = {delta_tot}")
    print("-" * 20)

    print("The problem asks for the thickness of 'the' double point, suggesting there is only one.")
    print(f"Therefore, the thickness of the double point is {delta_tot}.")

solve_curve_thickness()