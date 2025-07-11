def estimate_transition_temp(n, base_n=5, base_temp=25.0, slope=-4.0):
    """
    Estimates the transition temperature of a hypothetical CnH2n+1-Ph-CN liquid crystal
    based on the number of carbons (n) in the alkyl chain.

    This model is based on the design principles provided:
    - A base case of n=5 is assumed to give a transition temperature near room temp (~25 C).
    - As per the prompt, increasing chain length is assumed to decrease the temperature.
    - As per the prompt, decreasing chain length is assumed to increase the temperature.

    Args:
        n (int): The number of carbon atoms in the alkyl chain.
        base_n (int): The baseline number of carbons.
        base_temp (float): The baseline transition temperature in Celsius.
        slope (float): The change in temperature per unit change in n.

    Returns:
        float: The estimated transition temperature in Celsius.
    """
    # The equation to estimate temperature based on chain length n
    temperature = base_temp + (n - base_n) * slope
    return temperature

def main():
    """
    Main function to demonstrate the liquid crystal design principle.
    """
    print("Liquid Crystal Design: Estimating Transition Temperature")
    print("Based on the general structure: C(n)H(2n+1)-Ph-CN")
    print("Model based on design rule: 'start with a pentyl chain (n=5)' for ~25 C.")
    print("----------------------------------------------------------------------\n")

    # Case 1: Baseline from the design guide (n=5)
    n1 = 5
    temp1 = estimate_transition_temp(n1)
    print(f"For a starting alkyl chain with n = {n1} (pentyl):")
    # Output each number in the final equation
    print(f"Equation: T = 25.0 + ({n1} - 5) * (-4.0)")
    print(f"Estimated Transition Temperature: {temp1:.1f} C\n")

    # Case 2: Demonstrate Rule F: "If transition temperature is too high, increase chain length"
    # To demonstrate this, let's assume our target is 15 C and 25 C is too high.
    n2 = 8
    temp2 = estimate_transition_temp(n2)
    print(f"Demonstrating Rule F: If temperature is too high, increase chain length.")
    print(f"For an alkyl chain with n = {n2} (octyl):")
    # Output each number in the final equation
    print(f"Equation: T = 25.0 + ({n2} - 5) * (-4.0)")
    print(f"Estimated Transition Temperature: {temp2:.1f} C\n")

    # Case 3: Demonstrate Rule F: "If too low, decrease chain length"
    # To demonstrate this, let's assume our target is 35 C and 25 C is too low.
    n3 = 3
    temp3 = estimate_transition_temp(n3)
    print(f"Demonstrating Rule F: If temperature is too low, decrease chain length.")
    print(f"For an alkyl chain with n = {n3} (propyl):")
    # Output each number in the final equation
    print(f"Equation: T = 25.0 + ({n3} - 5) * (-4.0)")
    print(f"Estimated Transition Temperature: {temp3:.1f} C\n")

if __name__ == "__main__":
    main()