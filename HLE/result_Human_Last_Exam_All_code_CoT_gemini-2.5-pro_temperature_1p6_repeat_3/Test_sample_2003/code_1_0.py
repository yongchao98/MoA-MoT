def derive_and_show_sharps_formula():
    """
    This function demonstrates the derivation of the total number of sharps
    in the key signatures of 12 chromatically adjacent keys.

    The logic shows that the sum is constant regardless of the starting note,
    as the set of key signatures is always a complete chromatic collection.
    The number of sharps for each of the 12 unique major keys (when expressed
    only with sharps) corresponds to a permutation of the integers from 0 to 11.
    Therefore, the total sum is the sum of integers from 0 to 11.
    """

    print("The final formula for the sum of sharps (T) is derived by summing the number of sharps for all 12 unique major key signatures.")
    print("The number of sharps for these 12 keys are a permutation of the numbers 0 through 11.")
    print("\nThe final equation is the sum of these numbers:")
    
    # The numbers to sum are the integers from 0 to 11.
    numbers = list(range(12))
    
    # Calculate the sum.
    total = sum(numbers)
    
    # Create the full equation string.
    equation_str = " + ".join(map(str, numbers)) + f" = {total}"
    
    print(equation_str)
    
    print("\nThis simplifies to the formula:")
    print("T = 66")

# Execute the function to display the derivation.
derive_and_show_sharps_formula()