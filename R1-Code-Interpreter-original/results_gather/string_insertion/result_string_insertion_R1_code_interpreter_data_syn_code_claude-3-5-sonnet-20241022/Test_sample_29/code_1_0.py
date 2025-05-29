def apply_T_operation(s):
    # Define the patterns and their corresponding insertions
    patterns = {
        'ABCD': 'A',
        'BCDE': 'B',
        'CDEA': 'C',
        'DEAB': 'D',
        'EABC': 'E'
    }
    
    # Convert to list for easier manipulation
    result = list(s)
    insertions = []
    
    # Print each substring check for verification
    for i in range(len(s)-3):
        substr = s[i:i+4]
        print(f"Checking position {i}: {substr}")
        if substr in patterns:
            insertions.append((i+4, patterns[substr]))
            print(f"Found pattern {substr}, will insert {patterns[substr]}")
    
    # Apply insertions from left to right
    offset = 0
    for pos, char in sorted(insertions):
        result.insert(pos + offset, char)
        offset += 1
    
    return ''.join(result)

input_str = "CDEAAEABCCBCDE"
result = apply_T_operation(input_str)
print("\nFinal result:", result)