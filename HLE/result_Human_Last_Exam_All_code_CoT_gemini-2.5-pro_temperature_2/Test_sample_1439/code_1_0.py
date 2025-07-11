def find_critical_exponent_order():
    """
    This script explains and determines the order in the coupling constant 'u'
    at which the critical exponent ν for the φ⁴ theory acquires its first
    non-vanishing contribution.
    """
    print("This script will determine the order of the first correction to the critical exponent ν.")
    print("-" * 75)

    # 1. State the context and the mean-field (zeroth-order) result.
    print("In φ⁴ scalar field theory, the critical exponent ν is analyzed using a perturbative")
    print("expansion in the coupling constant 'u' around d=4 spacetime dimensions.")
    print("\nAt the zeroth order (in mean-field theory), the exponent has the value ν = 1/2.")
    print("This means 1/ν = 2.")

    # 2. Present the one-loop expansion formula for 1/ν.
    print("\nThe interactions introduce corrections. The one-loop calculation provides an expansion for 1/ν:")

    # Define the components of the equation to be printed.
    lhs_numerator = 1
    lhs_denominator = "ν"
    zeroth_order_term = 2
    first_order_coeff_numerator = "N+2"
    first_order_coeff_denominator = 6
    coupling_constant = "u"
    coupling_power = 1

    # Print the equation, referencing each number.
    print("\nThe one-loop formula is:")
    print(f"  {lhs_numerator}/{lhs_denominator} = {zeroth_order_term} - ({first_order_coeff_numerator}/{first_order_coeff_denominator}) * {coupling_constant}^{coupling_power} + O({coupling_constant}^2)")
    print("\n(Where N is the number of components of the scalar field.)")
    print("-" * 75)


    # 3. Analyze the formula and identify the first correction term.
    print("Let's analyze the terms in this equation:")
    print(f"* The first term on the right, '{zeroth_order_term}', is the constant zeroth-order term corresponding to the mean-field result.")
    print(f"* The second term, '-({first_order_coeff_numerator}/{first_order_coeff_denominator}) * {coupling_constant}^{coupling_power}', is the first correction term.")

    # 4. State the conclusion based on the power of 'u' in the first correction.
    print(f"\nThis first correction is proportional to the coupling constant '{coupling_constant}' raised to the power of {coupling_power}.")
    print("\nTherefore, the initial non-vanishing contribution to ν appears at the first order in 'u'.")

find_critical_exponent_order()
<<<1>>>