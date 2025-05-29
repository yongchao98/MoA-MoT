def apply_T_operation(s):
    # Store the patterns and their corresponding insertions
    patterns = {
        'ABCD': 'A',
        'BCDE': 'B',
        'CDEA': 'C',
        'DEAB': 'D',
        'EABC': 'E'
    }
    
    result = list(s)
    i = 0
    while i <= len(result) - 4:
        substr = ''.join(result[i:i+4])
        if substr in patterns:
            result.insert(i+4, patterns[substr])
            i += 1  # Move to next position after insertion
        else:
            i += 1
    
    final_string = ''.join(result)
    print(final_string)

# Test with the given string
s = "AEABCBEABCBEABCE"
apply_T_operation(s)