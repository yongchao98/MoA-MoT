def generate_sequence_s4(num_terms):
    """
    Generates and prints terms of a sequence based on the rule
    s[n] = s[s[n-1]] + s[s[n-2]].
    """
    s = {}
    if num_terms >= 1:
        s[1] = 1
    if num_terms >= 2:
        s[2] = 1
    
    # Print base cases
    print("s[1] = 1")
    if num_terms >= 2:
        print("s[2] = 1")

    for n in range(3, num_terms + 1):
        # The rule for s[n] depends on the previous two terms s[n-1] and s[n-2]
        # These values are used as indices to access earlier terms in the sequence.
        
        # Check if the required indices are within the bounds of the generated sequence
        if s[n-1] in s and s[n-2] in s:
            val_sn_minus_1 = s[n-1]
            val_sn_minus_2 = s[n-2]
            
            # The indices for the sum are the values from the previous steps
            index1 = val_sn_minus_1
            index2 = val_sn_minus_2

            s[n] = s[index1] + s[index2]
            
            # As requested, output the numbers in the final equation for the first few terms
            if n <= 8:
                 print(f"s[{n}] = s[s[{n-1}]] + s[s[{n-2}]] = s[{val_sn_minus_1}] + s[{val_sn_minus_2}] = {s[index1]} + {s[index2]} = {s[n]}")
            else:
                 print(f"s[{n}] = {s[n]}")

        else:
            print(f"Cannot calculate s[{n}] due to missing previous terms.")
            break

# Generate and print the first 42 terms of the sequence
generate_sequence_s4(42)
