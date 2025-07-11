def display_fixed_point_derivation():
    """
    Displays the derivation of the Wilson-Fisher fixed point coupling
    and prints its leading order expression.
    """
    # Define constants and symbols for the equations
    numerator_coeff = 16
    denominator_coeff = 3
    pi_symbol = "π"
    coupling_symbol = "u"
    fixed_point_coupling_symbol = "u*"
    epsilon_symbol = "ε"

    print("Derivation of the Wilson-Fisher Fixed Point in ϕ⁴ Theory:")
    print("-" * 60)
    
    # Step 1: State the beta function
    print(f"1. The one-loop beta function for the dimensionless coupling '{coupling_symbol}' in d = 4 - {epsilon_symbol} dimensions is given by:")
    print(f"   β({coupling_symbol}) = -{epsilon_symbol}{coupling_symbol} + (3 * {coupling_symbol}²)/(16 * {pi_symbol}²)\n")

    # Step 2: Apply the fixed point condition
    print(f"2. A fixed point, denoted by '{fixed_point_coupling_symbol}', occurs when the beta function is zero, i.e., β({fixed_point_coupling_symbol}) = 0.")
    print(f"   This gives the equation: 0 = -{epsilon_symbol}*{fixed_point_coupling_symbol} + (3 * {fixed_point_coupling_symbol}²)/(16 * {pi_symbol}²)\n")

    # Step 3: Solve for the non-trivial fixed point
    print(f"3. We seek the non-trivial solution where {fixed_point_coupling_symbol} ≠ 0. We can rearrange the equation:")
    print(f"   {epsilon_symbol}*{fixed_point_coupling_symbol} = (3 * {fixed_point_coupling_symbol}²)/(16 * {pi_symbol}²)")
    print(f"   Dividing both sides by {fixed_point_coupling_symbol}:")
    print(f"   {epsilon_symbol} = (3 * {fixed_point_coupling_symbol})/(16 * {pi_symbol}²)\n")

    # Step 4: Isolate the fixed point coupling to get the final expression
    print(f"4. Solving for '{fixed_point_coupling_symbol}' gives the leading order expression for the fixed point coupling:")
    print(f"\n   Final Expression: {fixed_point_coupling_symbol} = ({numerator_coeff} * {pi_symbol}² / {denominator_coeff}) * {epsilon_symbol}\n")

    # Highlight the numbers in the final equation as requested
    print("In this final equation:")
    print(f"- The number in the numerator is: {numerator_coeff}")
    print(f"- The term in the numerator also includes: {pi_symbol}²")
    print(f"- The number in the denominator is: {denominator_coeff}")
    print(f"- '{epsilon_symbol}' is the small parameter, defined as {epsilon_symbol} = 4 - d (where d is the number of spacetime dimensions).")

# Execute the function to display the result
if __name__ == "__main__":
    display_fixed_point_derivation()