def calculate_purification_product(F1, F2):
    """
    Calculates the product of the successful output fidelity and the success probability
    for the described GHZ state purification protocol.

    Args:
        F1 (float): The fidelity of the input 3-qubit GHZ state.
        F2 (float): The fidelity of the input 2-qubit Bell state.

    Returns:
        float: The product of the success probability and the output fidelity.
    """
    # The analytically derived formula for the product is: F1 * (2*F2 + 1) / 3
    # Check for valid fidelity ranges
    if not (0 <= F1 <= 1 and 0 <= F2 <= 1):
        raise ValueError("Fidelities F1 and F2 must be between 0 and 1.")

    # Numerator of the formula
    numerator = F1 * (2 * F2 + 1)
    
    # Denominator of the formula
    denominator = 3
    
    result = numerator / denominator
    
    # Print the final equation with the numbers substituted, as requested.
    print(f"The calculation for the product is based on the formula: F1 * (2 * F2 + 1) / 3")
    print(f"Substituting F1 = {F1} and F2 = {F2}:")
    print(f"Product = ({F1} * (2 * {F2} + 1)) / 3")
    print(f"Product = {result}")
    
    return result

if __name__ == '__main__':
    # You can change these values to your specific fidelities.
    # Example values are provided.
    input_fidelity_F1 = 0.95  # Fidelity of the input GHZ state
    input_fidelity_F2 = 0.9   # Fidelity of the input Bell state

    try:
        calculate_purification_product(input_fidelity_F1, input_fidelity_F2)
    except ValueError as e:
        print(f"Error: {e}")
