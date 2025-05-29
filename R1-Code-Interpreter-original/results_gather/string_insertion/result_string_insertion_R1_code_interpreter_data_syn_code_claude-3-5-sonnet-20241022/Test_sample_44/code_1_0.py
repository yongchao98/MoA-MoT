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
    i = 0
    
    # Keep track of how many characters we've inserted
    offset = 0
    
    # Scan the string
    while i <= len(s) - 4:
        # Get substring of length 4 starting at position i
        substr = s[i:i+4]
        if substr in patterns:
            # Insert the character at the correct position
            result.insert(i + 4 + offset, patterns[substr])
            offset += 1
        i += 1
    
    return ''.join(result)

# Test the string
input_str = "CCDEAACDEAEDEAB"
result = apply_T_operation(input_str)
print(result)