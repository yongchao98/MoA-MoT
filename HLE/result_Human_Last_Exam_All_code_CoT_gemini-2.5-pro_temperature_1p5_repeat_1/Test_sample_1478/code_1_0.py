import math

# This class simulates a term in a Titan expression, like '24/5e3'.
# It does not handle expression lists, only individual terms.
class TitanTerm:
    def __init__(self, num, den=1, exp=0):
        # A Titan register holds expressions. We check validity at each step.
        if not (0 <= num <= 63 and 0 < den <= 63):
            # This check is illustrative; the real check happens during operations.
            pass
        self.num = num
        self.den = den
        self.exp = exp

    def __repr__(self):
        return f"{self.num}/{self.den}e{self.exp}"

# This function simulates the multiplication of two Titan terms.
# It returns a new term or raises an error if the product exceeds the 6-bit limit.
def titan_multiply_terms(term1, term2, verbose=True):
    # This is the crucial check: the product of numerators cannot exceed 63.
    if term1.num * term2.num > 63:
        if verbose:
            print(f"OPERATION FAILED: Product of numerators {term1.num} * {term2.num} = {term1.num * term2.num} is > 63.")
        raise ValueError("Numerator out of 6-bit range")
    
    # The product of denominators also cannot exceed 63 without simplification.
    if term1.den * term2.den > 63:
        # Note: A real implementation would try to simplify this fraction.
        # For this problem, we show the constraint.
        if verbose:
            print(f"OPERATION FAILED: Product of denominators {term1.den} * {term2.den} = {term1.den * term2.den} is > 63.")
        raise ValueError("Denominator out of 6-bit range")

    new_num = term1.num * term2.num
    new_den = term1.den * term2.den
    new_exp = term1.exp + term2.exp
    
    # Here, we would simplify the resulting fraction (e.g., using GCD).
    # For simplicity, we assume we can handle valid products.
    return TitanTerm(new_num, new_den, new_exp)

# Main simulation logic
def simulate_calculation():
    print("Titan Feasibility Study: Calculating Gravitational Force")
    print("=" * 50)
    print("Objective: Calculate F = G * M * m / r^2")
    print("First, we must calculate the mass M = rho * 4/3 * pi * R^3")
    print("-" * 50)
    print("Step 1: Define constants as Titan terms")
    # rho = 1200 kg/m^3 = 1.2e3 = 6/5 e3
    rho = TitanTerm(6, 5, 3)
    # 4/3
    four_thirds = TitanTerm(4, 3)
    # pi approx 22/7, which we expand as (3/1 + 1/7) to manage intermediate products
    pi_whole = TitanTerm(3, 1)
    pi_frac = TitanTerm(1, 7)
    # R = 2e6 m, so R^3 = 8e18 m^3
    R_cubed = TitanTerm(8, 1, 18)
    
    print(f"rho = {rho}, 4/3 = {four_thirds}, pi approx = ({pi_whole} + {pi_frac}), R^3 = {R_cubed}\n")

    print("Step 2: Calculate M by evaluating the expression piece by piece")
    # Let's calculate C1 = 4/3 * rho
    print(f"Attempting C1 = {four_thirds} * {rho}")
    try:
        # (4/3) * (6/5) = 24/15. Simplify by GCD(24,15)=3 -> 8/5
        C1 = TitanTerm(8, 5, 3) 
        print(f"SUCCESS: C1 = 8/5e3. The expression is valid.\n")
    except ValueError as e:
        print(e)
        return False
        
    # Now, multiply C1 by pi. According to Titan rules, we must multiply by each part of pi's expression.
    # C1 * (pi_whole + pi_frac) = (C1 * pi_whole) + (C1 * pi_frac)
    # Part 1: C1 * pi_whole
    print(f"Attempting C2_part1 = C1 * pi_whole = {C1} * {pi_whole}")
    try:
        # 8/5e3 * 3/1 = 24/5e3. Both 8*3=24 and 5*1=5 are valid.
        C2_part1 = titan_multiply_terms(C1, pi_whole)
        print(f"SUCCESS: C2_part1 = {C2_part1}. The term is valid.")
    except ValueError:
        return False
    
    # Part 2: C1 * pi_frac
    print(f"Attempting C2_part2 = C1 * pi_frac = {C1} * {pi_frac}")
    try:
        # 8/5e3 * 1/7 = 8/35e3. Both 8*1=8 and 5*7=35 are valid.
        C2_part2 = titan_multiply_terms(C1, pi_frac)
        print(f"SUCCESS: C2_part2 = {C2_part2}. The term is valid.\n")
    except ValueError:
        return False

    print("The expression for Mass so far (without R^3) is represented by the sum:")
    print(f"M_intermediate = ({C2_part1}) + ({C2_part2})")
    print("-" * 50)
    
    print("Step 3: Multiply the intermediate mass expression by R_cubed")
    print("This requires multiplying each term in the expression by R_cubed.")
    print(f"Attempting: ({C2_part1}) * {R_cubed}")
    
    try:
        # The operation is (24/5 e3) * (8/1 e18).
        # It requires calculating the new numerator: 24 * 8.
        print("Checking numerator product: 24 * 8 = 192")
        titan_multiply_terms(C2_part1, R_cubed)
        # This line will not be reached because the function above raises an error.
    except ValueError:
        print("\nCONCLUSION: The calculation is NOT feasible.")
        print("The operation requires multiplying 24 by 8, resulting in 192.")
        print("192 cannot be represented by a 6-bit integer (0-63).")
        print("As the planet's mass cannot be calculated within Titan's constraints, the force cannot be determined.")
        return False

    return True

if not simulate_calculation():
    print("\nFinal Answer: The specified computation is not possible on the Titan architecture.")
    print('<<<N0>>>')
