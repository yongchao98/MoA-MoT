def apply_T(s):
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
    i = 0
    
    while i <= len(result) - 4:
        substring = ''.join(result[i:i+4])
        if substring in patterns:
            # Insert the corresponding character
            result.insert(i + 4, patterns[substring])
            i += 5  # Skip past the inserted character
        else:
            i += 1
    
    return ''.join(result)

# Test with the given string
input_string = "ABCDBABCDCCDEA"
result = apply_T(input_string)
print(result)