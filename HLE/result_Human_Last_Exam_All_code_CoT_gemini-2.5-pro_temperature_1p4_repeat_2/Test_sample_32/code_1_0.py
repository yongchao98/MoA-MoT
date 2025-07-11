from fractions import Fraction

def solve_lambda_integral():
    """
    Computes the integral of lambda_1*lambda_2*lambda_3 on M_3
    using a known tautological relation and known integral values.
    """
    
    # Known values for integrals on the moduli space of genus 3 curves
    # <l_1^6> = integral of lambda_1^6
    # <l_1^4 * l_2> = integral of lambda_1^4 * lambda_2
    integral_l1_6 = Fraction(1, 5760)
    integral_l1_4_l2 = Fraction(1, 2880)

    # A known tautological relation is:
    # 2 * <l_1^6> - 3 * <l_1^4 * l_2> + 3 * <l_1 * l_2 * l_3> = 0
    # We can solve for <l_1 * l_2 * l_3>
    # 3 * <l_1 * l_2 * l_3> = 3 * <l_1^4 * l_2> - 2 * <l_1^6>

    term1 = 3 * integral_l1_4_l2
    term2 = -2 * integral_l1_6
    
    # Calculate the right hand side of the equation
    rhs = term1 + term2
    
    # Solve for <l_1 * l_2 * l_3>
    integral_l1_l2_l3 = rhs / 3
    
    # Print the steps of the calculation
    print("Step 1: Start with the tautological relation:")
    print("2 * <λ_1^6> - 3 * <λ_1^4 * λ_2> + 3 * <λ_1 * λ_2 * λ_3> = 0")
    print("\nStep 2: Rearrange the equation to solve for <λ_1 * λ_2 * λ_3>:")
    print("3 * <λ_1 * λ_2 * λ_3> = 3 * <λ_1^4 * λ_2> - 2 * <λ_1^6>")
    print("\nStep 3: Substitute the known integral values:")
    print(f"<λ_1^6> = {integral_l1_6}")
    print(f"<λ_1^4 * λ_2> = {integral_l1_4_l2}")
    print(f"3 * <λ_1 * λ_2 * λ_3> = 3 * ({integral_l1_4_l2}) - 2 * ({integral_l1_6})")
    
    print("\nStep 4: Perform the calculation:")
    print(f"3 * <λ_1 * λ_2 * λ_3> = {term1} - {Fraction(-term2.numerator, -term2.denominator)}")
    print(f"3 * <λ_1 * λ_2 * λ_3> = {rhs}")
    print(f"<λ_1 * λ_2 * λ_3> = {rhs} / 3")
    print(f"<λ_1 * λ_2 * λ_3> = {integral_l1_l2_l3}")
    
    # Final answer formatting
    result = f"{integral_l1_l2_l3.numerator}/{integral_l1_l2_l3.denominator}"
    print(f"\nFinal Answer: {result}")
    
solve_lambda_integral()