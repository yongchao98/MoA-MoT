def apply_T_operation(s):
    # Define the patterns and their corresponding insertions
    patterns = {
        'ABCD': 'A',
        'BCDE': 'B',
        'CDEA': 'C',
        'DEAB': 'D',
        'EABC': 'E'
    }
    
    result = list(s)
    insertions = []
    
    # Scan the string
    i = 0
    while i <= len(s) - 4:
        substr = s[i:i+4]
        if substr in patterns:
            insertions.append((i+4, patterns[substr]))
            i += 4  # Move past the current pattern
        else:
            i += 1
    
    # Apply insertions from left to right
    offset = 0
    for pos, char in sorted(insertions):
        result.insert(pos + offset, char)
        offset += 1
    
    return ''.join(result)

# Test the function with the given string
input_string = "ABCDDCDEAEDEABD"
result = apply_T_operation(input_string)
print(result)