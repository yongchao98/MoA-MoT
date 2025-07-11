import collections

def solve_cellular_automaton():
    """
    This function deduces the possible elementary cellular automaton rules that could
    produce the pattern shown in the image.
    """
    # Step 1: Analyze the visual pattern to deduce constraints on the rule.
    # The image shows a checkerboard pattern growing from a single central cell.
    # By observing which 3-cell neighborhoods produce which output cell in the next
    # generation, we can determine parts of the rule.
    # A '1' represents a black cell and a '0' represents a white cell.
    #
    # Example transitions observed:
    # - At the edge of the pattern, a neighborhood of (0,0,1) results in a 1.
    # - In a sea of white cells, (0,0,0) results in a 0.
    # - A single '1' surrounded by '0's, i.e., (0,1,0), results in a 0.
    # - In the checkerboard area, a neighborhood of (1,0,1) results in a 1.
    #
    # From this visual analysis, we find the following required mappings:
    constraints = {
        # neighborhood: resulting_state
        (0, 0, 0): 0,
        (0, 0, 1): 1,
        (0, 1, 0): 0,
        (1, 0, 0): 1,
        (1, 0, 1): 1,
    }

    # The neighborhoods (1,1,1), (1,1,0), and (0,1,1) do not appear in this specific
    # evolution, so their output is not constrained by the image.
    
    print("This program finds all elementary cellular automaton rules consistent with the provided image.")
    print("The analysis reveals the following required transitions:")
    for neighborhood, result in constraints.items():
        print(f"Neighborhood {neighborhood} must result in {result}")

    # Step 2: Construct all possible rules based on the constraints.
    # An elementary rule is defined by 8 bits, for the 8 possible neighborhoods.
    # The Wolfram convention numbers these neighborhoods from '111' down to '000'.
    # The rule's integer value is the decimal of this 8-bit binary string.
    #
    # Bits are numbered 7 down to 0:
    # Bit 7: (1,1,1) - Unconstrained
    # Bit 6: (1,1,0) - Unconstrained
    # Bit 5: (1,0,1) -> 1
    # Bit 4: (1,0,0) -> 1
    # Bit 3: (0,1,1) - Unconstrained
    # Bit 2: (0,1,0) -> 0
    # Bit 1: (0,0,1) -> 1
    # Bit 0: (0,0,0) -> 0
    
    # Calculate the base value from the fixed (constrained) bits.
    # The rule's binary form is: ? ? 1 1 ? 0 1 0
    # Base value = 00110010 (binary)
    base_value = (1 << 5) + (1 << 4) + (0 << 2) + (1 << 1) + (0 << 0)
    
    print(f"\nThe determined part of the rule is '..11.010', which corresponds to a base value of {base_value}.")
    
    # The unconstrained bits correspond to powers of 2.
    unconstrained_bit_values = [1 << 7, 1 << 6, 1 << 3]
    
    print(f"The unconstrained bits are at positions 7, 6, and 3, with values {unconstrained_bit_values[0]}, {unconstrained_bit_values[1]}, and {unconstrained_bit_values[2]}.")
    print("We must find all combinations of these values added to the base value.")
    
    possible_rules = []
    # There are 2^3 = 8 combinations for the 3 unconstrained bits.
    for i in range(8):
        rule = base_value
        # Use the binary representation of i to decide whether to add each unconstrained value.
        # Example calculation for i=7 (binary 111), which yields the largest rule number.
        if i == 7:
            print("\nExample calculation for the largest rule number (all unconstrained bits set to 1):")
            calculation_str = f"{base_value}"
        
        if (i >> 2) & 1:  # Corresponds to bit 7
            rule += unconstrained_bit_values[0]
            if i == 7: calculation_str += f" + {unconstrained_bit_values[0]}"
        if (i >> 1) & 1:  # Corresponds to bit 6
            rule += unconstrained_bit_values[1]
            if i == 7: calculation_str += f" + {unconstrained_bit_values[1]}"
        if (i >> 0) & 1:  # Corresponds to bit 3
            rule += unconstrained_bit_values[2]
            if i == 7: calculation_str += f" + {unconstrained_bit_values[2]}"
            
        if i == 7:
            print(f"{calculation_str} = {rule}")
            
        possible_rules.append(rule)

    possible_rules.sort()
    
    print("\nThe complete list of possible rules, sorted in increasing order, is:")
    print(','.join(map(str, possible_rules)))

solve_cellular_automaton()