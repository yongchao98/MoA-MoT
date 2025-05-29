def apply_T_operation(s):
    # Store the patterns and their corresponding insertions
    patterns = {
        'ABCD': 'A',
        'BCDE': 'B',
        'CDEA': 'C',
        'DEAB': 'D',
        'EABC': 'E'
    }
    
    # Find all occurrences of all patterns
    insertions = []
    for i in range(len(s)):
        substr = s[i:i+4]
        if substr in patterns:
            print(f"Found pattern {substr} at position {i}, will insert {patterns[substr]}")
        if substr in patterns:
            insertions.append((i+4, patterns[substr]))
    
    # Sort insertions by position
    insertions.sort()
    print("\nInsertions to be made:", insertions)
    
    # Apply insertions from left to right
    result = list(s)
    offset = 0
    for pos, char in insertions:
        result.insert(pos + offset, char)
        offset += 1
    
    final_string = ''.join(result)
    print("\nFinal string:", final_string)

# Test with the given string
s = "BCDEEEABCCEABCC"
apply_T_operation(s)