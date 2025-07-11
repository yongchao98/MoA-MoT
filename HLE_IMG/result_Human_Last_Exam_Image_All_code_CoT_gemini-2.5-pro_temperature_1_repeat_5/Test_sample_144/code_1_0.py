def solve_mimicry_puzzle():
    """
    Solves the insect mimicry puzzle by identifying and pairing mimics with damage-causers.
    """
    
    # The pairs are determined by matching the appearance of the mimic
    # with the type of leaf damage caused by the causer.
    # Mimic A (adult beetle) mimics the linear damage of its larva (Causer B).
    # Mimic C (moth) is paired with Causer D (adult beetle) by elimination.
    # Mimic E (leaf insect) mimics the chewing damage of the Causer F (katydid).
    
    pairs = [('A', 'B'), ('C', 'D'), ('E', 'F')]
    
    # Format the pairs into the required string format "AB, CD, EF".
    formatted_pairs = []
    for mimic, causer in pairs:
        formatted_pairs.append(f"{mimic}{causer}")
        
    result_string = ", ".join(formatted_pairs)
    
    print("The matched pairs of mimic and damage-causing insects are:")
    print(result_string)

solve_mimicry_puzzle()