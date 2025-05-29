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
    
    # Print each substring of length 4 and whether it matches
    for i in range(len(s)-3):
        substr = s[i:i+4]
        if substr in patterns:
            print(f"Match found at position {i}: {substr} -> insert {patterns[substr]}")
            insertions.append((i+4, patterns[substr]))
    
    # Apply insertions from left to right
    offset = 0
    for pos, char in sorted(insertions):
        result.insert(pos + offset, char)
        offset += 1
        
    final = ''.join(result)
    print(f"\nInsertions (position, char): {insertions}")
    return final

s = "DEABBEABCCABCDA"
result = apply_T_operation(s)
print(f"\nFinal result: {result}")