def apply_T_operation(s):
    # Store all patterns and their corresponding insertions
    patterns = {
        'ABCD': 'A',
        'BCDE': 'B',
        'CDEA': 'C',
        'DEAB': 'D',
        'EABC': 'E'
    }
    
    # Convert string to list for easier manipulation
    result = list(s)
    i = 0
    insertions = []
    
    # Scan string for patterns
    while i <= len(s) - 4:
        substr = s[i:i+4]
        if substr in patterns:
            insertions.append((i+4, patterns[substr]))
        i += 1
    
    # Apply insertions from left to right
    offset = 0
    for pos, char in insertions:
        result.insert(pos + offset, char)
        offset += 1
    
    final_string = ''.join(result)
    print(final_string)

# Test with the given string
s = "ACEABCBEABCBBCDEE"
apply_T_operation(s)