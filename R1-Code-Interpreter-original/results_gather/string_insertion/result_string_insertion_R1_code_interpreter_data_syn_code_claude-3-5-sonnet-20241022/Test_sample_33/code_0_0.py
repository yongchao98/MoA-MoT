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
    insertions = 0  # Track number of insertions to adjust index
    
    while i < len(s) - 3:  # We need at least 4 characters to match
        # Get the current 4-character substring
        substr = s[i:i+4]
        
        # Check if it matches any pattern
        if substr in patterns:
            # Insert the corresponding character at the correct position
            result.insert(i + 4 + insertions, patterns[substr])
            insertions += 1
        
        i += 1
    
    # Convert back to string
    result = ''.join(result)
    print(result)

# Test with the given string
s = "DEABEEABCADEAB"
apply_T_operation(s)