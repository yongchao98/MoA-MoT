def apply_T_operation(s):
    # Define patterns and their corresponding insertions
    patterns = {
        'ABCD': 'A',
        'BCDE': 'B',
        'CDEA': 'C',
        'DEAB': 'D',
        'EABC': 'E'
    }
    
    # Convert to list for easier manipulation
    result = list(s)
    offset = 0  # Track position shifts due to insertions
    
    # Scan through original string positions
    for i in range(len(s)-3):  # -3 because patterns are 4 chars long
        # Get current 4-char substring from original position
        substr = s[i:i+4]
        
        # Check if substring matches any pattern
        if substr in patterns:
            # Insert corresponding character at position after pattern
            result.insert(i + 4 + offset, patterns[substr])
            offset += 1
    
    # Convert back to string and return
    print(''.join(result))

# Test with given string
apply_T_operation("CADEABACDEADCDEAED")