import numpy as np

def calculate_potential():
    """
    Calculates the electric potential based on the derived formula and prints the results.
    """
    # Define the parameters for the calculation
    sigma_0 = 1e-9  # Surface charge density amplitude (C/m^2)
    k = 100       # Wave number (rad/m)
    a = 0.01      # Distance to upper plate (m)
    b = 0.01      # Distance to lower plate (m)
    epsilon_1 = 2.2 * 8.854e-12 # Permittivity of region 1 (F/m)
    epsilon_2 = 4.4 * 8.854e-12 # Permittivity of region 2 (F/m)

    # Point (x, y) where the potential is to be calculated
    x = 0.005 # m
    y = 0.005 # m (This is in the region 0 < y < a)

    print("The electric potential Φ(x, y) is determined by solving Laplace's equation with the given boundary conditions.")
    print("The correct derived solution corresponds to choice A.")
    
    print("\nThe full potential Φ(x,y) is:")
    print("Φ(x, y) =")
    print("  for 0 < y < a: [-σ₀ sinh(kb) sinh(k(y-a)) sin(kx)] / [k(ε₂cosh(ka)sinh(kb) + ε₁sinh(ka)cosh(kb))]")
    print("  for -b < y < 0: [ σ₀ sinh(ka) sinh(k(y+b)) sin(kx)] / [k(ε₂cosh(ka)sinh(kb) + ε₁sinh(ka)cosh(kb))]")
    
    print("\n--- Numerical Calculation ---")
    print(f"Given parameters:")
    print(f"σ₀ (sigma_0) = {sigma_0:.2e} C/m^2")
    print(f"k = {k} rad/m")
    print(f"a = {a} m")
    print(f"b = {b} m")
    print(f"ε₁ (epsilon_1) = {epsilon_1:.2e} F/m")
    print(f"ε₂ (epsilon_2) = {epsilon_2:.2e} F/m")
    print(f"Point of interest (x, y) = ({x}, {y}) m")

    if 0 < y < a:
        # Helper functions
        sinh = np.sinh
        cosh = np.cosh
        sin = np.sin
        
        # Calculate intermediate terms
        skb = sinh(k * b)
        sk_ya = sinh(k * (y - a))
        skx = sin(k * x)
        cka = cosh(k * a)
        ska = sinh(k * a)
        ckb = cosh(k * b)

        # Numerator
        numerator = -sigma_0 * skb * sk_ya * skx
        # Denominator
        denominator = k * (epsilon_2 * cka * skb + epsilon_1 * ska * ckb)
        # Final Potential
        potential = numerator / denominator

        print("\nCalculating the potential in the region 0 < y < a:")
        print("Final Equation with numbers:")
        print(f"Φ({x}, {y}) = [ -({sigma_0:.2e}) * sinh({k}*{b}) * sinh({k}*({y} - {a})) * sin({k}*{x}) ] / "
              f"\\\n             [ {k} * (({epsilon_2:.2e})*cosh({k}*{a})*sinh({k}*{b}) + ({epsilon_1:.2e})*sinh({k}*{a})*cosh({k}*{b})) ]")
        
        print(f"\nThe calculated potential at ({x}, {y}) is: {potential:.4f} V")
        
    elif -b < y < 0:
        # Helper functions
        sinh = np.sinh
        cosh = np.cosh
        sin = np.sin
        
        # Calculate intermediate terms
        ska = sinh(k * a)
        sk_yb = sinh(k * (y + b))
        skx = sin(k * x)
        cka = cosh(k * a)
        skb = sinh(k * b)
        ckb = cosh(k * b)

        # Numerator
        numerator = sigma_0 * ska * sk_yb * skx
        # Denominator
        denominator = k * (epsilon_2 * cka * skb + epsilon_1 * ska * ckb)
        # Final Potential
        potential = numerator / denominator
        
        print("\nCalculating the potential in the region -b < y < 0:")
        print("Final Equation with numbers:")
        print(f"Φ({x}, {y}) = [ ({sigma_0:.2e}) * sinh({k}*{a}) * sinh({k}*({y} + {b})) * sin({k}*{x}) ] / "
              f"\\n             [ {k} * (({epsilon_2:.2e})*cosh({k}*{a})*sinh({k}*{b}) + ({epsilon_1:.2e})*sinh({k}*{a})*cosh({k}*{b})) ]")

        print(f"\nThe calculated potential at ({x}, {y}) is: {potential:.4f} V")
        
    else:
        print("\nThe point (x, y) is on a boundary conductor where potential is 0, or outside the defined regions.")


if __name__ == '__main__':
    calculate_potential()
