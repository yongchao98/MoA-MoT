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
    offset = 0  # Track position shifts due to insertions
    
    # Scan through original string positions
    for i in range(len(s)):
        # Get the current 4-character substring starting at position i
        if i + 4 <= len(s):
            substr = s[i:i+4]
            # Check if this substring matches any pattern
            if substr in patterns:
                # Insert the corresponding character at the correct position
                result.insert(i + 4 + offset, patterns[substr])
                offset += 1
    
    # Convert back to string and return
    final_result = ''.join(result)
    print(final_result)

# Test with the given string
s = "EBCDEEEABCCBCDE"
apply_T_operation(s)