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
    offset = 0  # Track how many characters we've inserted
    
    while i < len(s):
        for pattern, insert_char in patterns.items():
            # Check if pattern exists starting at position i
            if i + len(pattern) <= len(s):
                substring = s[i:i + len(pattern)]
                if substring == pattern:
                    # Insert the character at the correct position
                    insert_pos = i + len(pattern) + offset
                    result.insert(insert_pos, insert_char)
                    offset += 1
        i += 1
    
    # Convert back to string
    result = ''.join(result)
    print(result)

# Test with the given string
s = "DAEABCAEABCCEABC"
apply_T_operation(s)