import itertools

def solve():
    """
    Deduces the ECA rules that could generate the initial part of the given pattern.
    """
    # The pattern for the first 3 rows (t=0, 1, 2).
    # White=0, Black=1. Width is chosen to see the whole evolution cone.
    pattern = [
        [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],  # t=0
        [0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0],  # t=1
        [0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0],  # t=2
    ]

    # Rule bits are for neighborhoods 111, 110, 101, 100, 011, 010, 001, 000
    # Initialize rule with None to indicate 'unconstrained'
    rule = {f"{i:03b}": None for i in range(8)}

    # We only need to check transitions from t=0 to t=1 and t=1 to t=2
    # This covers the generation of the first three rows of the pattern.
    for t in range(2):
        parent_row = pattern[t]
        child_row = pattern[t+1]
        for i in range(1, len(parent_row) - 1):
            neighborhood = tuple(parent_row[i-1 : i+2])
            neighborhood_str = "".join(map(str, neighborhood))
            output = child_row[i]
            
            # If we see a neighborhood for the first time, record its output
            if rule[neighborhood_str] is None:
                rule[neighborhood_str] = output
            # If we see it again, check for consistency (although we know it's consistent here)
            elif rule[neighborhood_str] != output:
                # This part is not triggered by the first 3 rows
                print(f"Error: Contradiction found for neighborhood {neighborhood_str}")
                return

    # Find unconstrained bits
    unconstrained_bits = []
    constrained_bits = {}
    
    # Wolfram convention orders neighborhoods from 111 down to 000
    neighborhood_order = [f"{i:03b}" for i in range(7, -1, -1)]

    for neighborhood_str in neighborhood_order:
        bit_value = rule[neighborhood_str]
        if bit_value is None:
            unconstrained_bits.append(neighborhood_str)
        else:
            constrained_bits[neighborhood_str] = bit_value

    # Generate all possible rule binaries
    possible_rules_binary = []
    num_unconstrained = len(unconstrained_bits)
    
    # Iterate through all possibilities (0/1) for the unconstrained bits
    for p in itertools.product([0, 1], repeat=num_unconstrained):
        temp_rule = constrained_bits.copy()
        for i, bit_key in enumerate(unconstrained_bits):
            temp_rule[bit_key] = p[i]
        
        # Build the full 8-bit rule string
        binary_string = ""
        for n_str in neighborhood_order:
            binary_string += str(temp_rule[n_str])
        possible_rules_binary.append(binary_string)

    # Convert binary strings to integer rule numbers and sort them
    possible_rules_int = sorted([int(b, 2) for b in possible_rules_binary])

    # Print the result
    print("The following transitions are observed from the first 3 rows:")
    for n_str in neighborhood_order:
        if n_str not in unconstrained_bits:
            print(f"Neighborhood {n_str} -> {constrained_bits[n_str]}")

    print("\nThe following neighborhood(s) are not observed in the first 2 parent rows:")
    for n_str in unconstrained_bits:
        print(f"Neighborhood {n_str} -> ? (unconstrained)")
    
    print("\nThis leads to the following possible rules (integer representation):")
    print(','.join(map(str, possible_rules_int)))

solve()