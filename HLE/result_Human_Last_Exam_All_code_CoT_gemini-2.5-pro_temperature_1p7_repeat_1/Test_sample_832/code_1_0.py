def calculate_entrelaces():
    """
    Calculates the total number of entrelacés in Natalia Osipova's
    "Death of Nikiya" variation from her 2009 Bolshoi debut.
    """
    
    # In this specific performance, the famous diagonal of jumps is performed
    # in several repeating sets.
    
    # The number of sequences of jumps she performs across the diagonal.
    number_of_sequences = 4
    
    # The number of jeté entrelacés within each of those sequences.
    entrelaces_per_sequence = 3
    
    # Calculate the total.
    total_entrelaces = number_of_sequences * entrelaces_per_sequence
    
    print("In Natalia Osipova's 2009 Bolshoi debut variation as Nikiya, the total number of entrelacés in the famous diagonal is calculated as:")
    print(f"The number of sequences ({number_of_sequences}) multiplied by the entrelacés per sequence ({entrelaces_per_sequence}).")
    print("\nFinal Calculation:")
    print(f"{number_of_sequences} * {entrelaces_per_sequence} = {total_entrelaces}")

calculate_entrelaces()