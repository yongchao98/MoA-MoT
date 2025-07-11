import math

# Helper class to represent fractions and perform calculations under Titan's rules.
class TitanFraction:
    def __init__(self, n, d=1):
        if not (0 <= n <= 31 and 0 <= d <= 31):
            # This check is for initial setup. The main constraint is on multiplication products.
            raise ValueError(f"Numerator {n} or Denominator {d} out of 5-bit range (0-31)")
        self.n = n
        self.d = d

    def __repr__(self):
        return f"{self.n}/{self.d}"

    # Custom multiplication that respects Titan's constraints
    def multiply(self, other, simplify_first=None):
        n1, d1 = self.n, self.d
        n2, d2 = other.n, other.d
        
        # Optional simplification before multiplication (e.g., (a/b)*(c/a) -> c/b)
        if simplify_first:
            common_factor_val = simplify_first[0]
            if n1 == common_factor_val and d2 == common_factor_val:
                n1 = 1
                d2 = 1
            elif n2 == common_factor_val and d1 == common_factor_val:
                n2 = 1
                d1 = 1
        
        # Check if the direct product exceeds the 5-bit limit
        if n1 * n2 > 31 or d1 * d2 > 31:
            raise ValueError(f"Operation ({n1}/{d1}) * ({n2}/{d2}) is illegal: "
                             f"Product {n1*n2} or {d1*d2} exceeds 31.")

        # If we are here, the operation is valid. Return the new unsimplified fraction.
        return TitanFraction(n1 * n2, d1 * d2)

    # Division is multiplication by the reciprocal
    def divide(self, other):
        reciprocal = TitanFraction(other.d, other.n)
        return self.multiply(reciprocal)
        
    def simplify(self):
        """Simplifies the fraction."""
        common_divisor = math.gcd(self.n, self.d)
        self.n //= common_divisor
        self.d //= common_divisor
        return self

def calculate_force_on_titan():
    """
    Simulates the calculation of the force F on the Titan computer.
    """
    print("--- Titan Computer Calculation ---")
    print("Goal: Calculate Force F = (d * m * g) / h\n")

    # 1. Define constants from the problem
    d = TitanFraction(20, 1)  # distance = 20m, aiming for center of lion
    h = TitanFraction(10, 1)  # height = 10m
    
    # Constants for mass calculation
    pi = TitanFraction(22, 7)    # π ≈ 22/7
    r = TitanFraction(1, 2)     # radius = 0.5cm = 1/2 cm
    rho = TitanFraction(9, 10)  # density = 0.9 kg/cm³ = 9/10 kg/cm³
    four_thirds = TitanFraction(4, 3)

    print(f"Initial Constants:")
    print(f"d = {d}, h = {h}, π = {pi}, r = {r}, ρ = {rho}\n")

    # 2. Calculate the rock's volume V = (4/3) * π * r^3 in cm³
    print("Step 1: Calculate Volume (V = 4/3 * π * r³)")
    r_cubed = TitanFraction(1, 8) # (1/2)³ = 1/8
    # V = (4/3) * (22/7) * (1/8). Order of operations matters.
    # (4/3)*(22/7) -> 88/21. Illegal (88 > 31).
    # We must reorder: (22/7)*(1/8) first.
    V_intermediate = pi.multiply(r_cubed) # (22/7)*(1/8) = 22/56. Illegal (56>31).
    
    # We must simplify before multiplying. Let's show a valid path.
    # (22/7) * (1/8) can be simplified as (11/7) * (1/4) = 11/28.
    print(f"To calculate V, we reorder and simplify: V = ((22/7) * (1/8)) * (4/3)")
    V_interim = TitanFraction(11, 28) # Manually simplified from 22/56
    print(f"  (22/7) * (1/8) -> simplified to {V_interim}")
    # Now V = (11/28) * (4/3). This is (11*4)/(28*3) = 44/84. Illegal.
    # We simplify again: (11/28)*(4/3) -> 11 * (4/28) / 3 -> 11 * (1/7) / 3 = 11/21
    V = TitanFraction(11, 21) # Manually simplified
    print(f"  ({V_interim}) * (4/3) -> simplified to V = {V}\n")

    # 3. Calculate mass m = V * ρ.
    print("Step 2: Calculate Mass (m = V * ρ)")
    print(f"  Attempting m = {V} * {rho} = ({V.n}*{rho.n})/({V.d}*{rho.d})")
    print(f"  This is illegal, as {V.n}*{rho.n} = {V.n*rho.n} which is > 31.")
    
    # Sacrifice precision to maintain constraint, as per rule #4.
    # V = 11/21 ≈ 0.52. Let's approximate it to 1/2.
    V_approx = TitanFraction(1, 2)
    print(f"  Approximation needed: V = {V} is approximated to {V_approx}")
    m = V_approx.multiply(rho) # (1/2)*(9/10) = 9/20. This is a valid operation.
    print(f"  m = {V_approx} * {rho} = {m} kg\n")

    # 4. Select g and calculate the final Force F.
    print("Step 3: Calculate Force (F = (d/h) * m * g)")
    # g ≈ 9.8 m/s². Approximations like 29/3 or 19/2 lead to illegal operations.
    # e.g., F = (9/10) * (29/3) -> (3/10)*(29/1) -> num = 87 > 31.
    # We must choose g = 10/1 for the calculation to be possible.
    g = TitanFraction(10, 1)
    print(f"  Chosen g = {g} (approx. 9.8) to allow calculation to proceed.\n")

    # 5. Execute the final calculation step-by-step
    print("Final Calculation Steps:")
    d_over_h = d.divide(h).simplify()
    print(f"  d/h = {d} / {h} = {d_over_h}")

    term2 = d_over_h.multiply(m, simplify_first=[2])
    # Manually simplify (2/1)*(9/20) -> (1/1)*(9/10)
    term2.simplify()
    print(f"  (d/h) * m = {d_over_h} * {m} -> simplified to {term2}")
    
    # Final multiplication: (9/10) * (10/1)
    # We simplify before multiplying the numerators and denominators.
    # 10/10 cancels out, leaving 9/1.
    final_F = term2.multiply(g, simplify_first=[10]).simplify()
    print(f"  ((d/h) * m) * g = {term2} * {g} -> simplified to {final_F}\n")

    print("--- Final Result ---")
    print("The final equation with the chosen fractional numbers is:")
    # Print the full equation
    print(f"F = (({d.n}/{d.d}) / ({h.n}/{h.d})) * ({V_approx.n}/{V_approx.d} * {rho.n}/{rho.d}) * ({g.n}/{g.d})")
    # Print the intermediate values
    print(f"F = ({d_over_h.n}/{d_over_h.d}) * ({m.n}/{m.d}) * ({g.n}/{g.d})")
    print(f"F = {final_F.n}/{final_F.d} Newtons")

calculate_force_on_titan()