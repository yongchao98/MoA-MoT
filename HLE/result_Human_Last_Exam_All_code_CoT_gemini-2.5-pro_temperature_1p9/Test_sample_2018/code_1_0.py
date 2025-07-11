def find_original_factor():
    """
    This function calculates and displays the value of the computational factor
    from the prior published work mentioned in the problem description.
    """
    
    # The computational factor in the original paper is given by the product
    # of a base value and a power of 10.
    
    base = 1.6
    exponent = 3
    
    # Calculate the final value of the computational factor
    value = base * (10**exponent)
    
    print(f"The equation for the computational factor is:")
    # We output each number used in the calculation as requested.
    print(f"{base} * 10^{exponent}")
    print(f"The value of the computational factor used in the prior work was: {value}")

find_original_factor()