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
    offset = 0  # Track how many characters we've inserted
    
    # Check each position in the original string
    for i in range(len(s)):
        # Get the substring starting at current position
        substr = s[i:i+4]
        
        # Check if any pattern matches
        if substr in patterns:
            # Insert the corresponding character at the correct position
            result.insert(i + 4 + offset, patterns[substr])
            offset += 1
    
    # Convert back to string
    final_result = ''.join(result)
    print(final_result)

# Test with the given string
s = "AAABCDBBCDEECDEAD"
apply_T_operation(s)