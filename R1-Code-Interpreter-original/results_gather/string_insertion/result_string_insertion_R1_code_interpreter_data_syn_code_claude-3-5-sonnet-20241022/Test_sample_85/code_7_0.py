def apply_T(s):
    # Initialize result string
    result = s[0]  # Start with first character
    
    # Process each character of original string
    for i in range(1, len(s)):
        # Add current character
        result += s[i]
        
        # Only check patterns when we have enough characters
        # and only check the original string positions
        if i >= 3:  # We have at least 4 characters to check
            original_pos = i - 3  # Start position in original string
            if original_pos + 4 <= len(s):  # Ensure we're checking original string pattern
                pattern = s[original_pos:original_pos + 4]
                if pattern == 'ABCD':
                    result += 'A'
                elif pattern == 'BCDE':
                    result += 'B'
                elif pattern == 'CDEA':
                    result += 'C'
                elif pattern == 'DEAB':
                    result += 'D'
                elif pattern == 'EABC':
                    result += 'E'
    
    print(result)

s = "EABCBDEABADEAB"
apply_T(s)