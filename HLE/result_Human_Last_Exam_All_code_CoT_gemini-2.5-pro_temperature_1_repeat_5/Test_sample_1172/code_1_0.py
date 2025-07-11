import sympy
from sympy import Symbol, Eq, pi

def calculate_inductance_change():
    """
    This script symbolically calculates the change in mutual inductance (ΔM)
    between two circuits when surrounded by magnetic concentrators.

    The final formula is constructed based on the physical derivation explained
    in the plan.
    """

    # Define the physical constants and geometric variables as symbols
    mu_0 = Symbol('μ₀', real=True, positive=True)  # Permeability of free space
    h = Symbol('h', real=True, positive=True)        # Wire separation in a circuit
    d = Symbol('d', real=True, positive=True)        # Distance between circuits
    R1 = Symbol('R₁', real=True, positive=True)      # Inner radius of the concentrator
    R2 = Symbol('R₂', real=True, positive=True)      # Outer radius of the concentrator
    delta_M = Symbol('ΔM')                           # The change in mutual inductance

    # Step 1: Expression for M₁, the mutual inductance for bare circuits
    # In the limit d >> h, M₁ ≈ (μ₀ * h²) / (2 * π * d²)
    M1 = (mu_0 * h**2) / (2 * pi * d**2)

    # Step 2 & 3: Expression for M₂, the mutual inductance with concentrators
    # M₂ = M₁ * (R₂ / R₁)
    M2 = M1 * (R2 / R1)

    # Step 4: Calculate the change, ΔM = M₂ - M₁
    final_expression_val = sympy.simplify(M2 - M1)

    # As requested, output the final equation showing each term/symbol
    print("The derived expression for the change in mutual inductance per unit length is:")
    # The equation is ΔM = (μ₀ * h²) / (2 * π * d²) * (R₂/R₁ - 1)
    # The following print statements build this equation piece by piece.
    
    # Numerator of the first term
    numerator = mu_0 * h**2
    # Denominator of the first term
    denominator = 2 * pi * d**2
    # Second term (the scaling factor from the concentrator)
    factor = (R2 / R1 - 1)

    print(f"ΔM = ({numerator}) / ({denominator}) * ({factor})")
    
    # Return the final symbolic expression for the answer block.
    return final_expression_val

if __name__ == '__main__':
    # Execute the calculation
    final_answer = calculate_inductance_change()
    
    # Format the final answer as requested
    pretty_answer = sympy.pretty(final_answer, use_unicode=True)
    print(f"\n<<<ΔM = {pretty_answer}>>>")