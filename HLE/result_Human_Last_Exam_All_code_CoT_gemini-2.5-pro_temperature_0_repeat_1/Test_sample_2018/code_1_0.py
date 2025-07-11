def find_computational_factor():
    """
    This function calculates and displays the value of the computational factor
    used in the prior published simulation-only work for the enthalpy-porosity method.
    The value is C = 1.6 * 10^6, as stated in the 1987 paper by Voller and Prakash.
    """
    
    # The components of the value from the scientific paper
    base_value = 1.6
    power_of = 10
    exponent = 6
    
    # Calculate the final value
    final_value = base_value * (power_of ** exponent)
    
    # Print the equation with each number
    print(f"The equation for the computational factor (C) is:")
    print(f"C = {base_value} * {power_of}^{exponent}")
    
    # Print the final calculated value
    print(f"\nThe calculated value of C is: {final_value}")

# Execute the function
find_computational_factor()