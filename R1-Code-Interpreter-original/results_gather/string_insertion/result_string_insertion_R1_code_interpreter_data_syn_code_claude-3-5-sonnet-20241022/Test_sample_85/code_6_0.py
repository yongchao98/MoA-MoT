def apply_T(s):
    # Process one character at a time and build result
    result = ""
    i = 0
    while i < len(s):
        # Add current character
        result += s[i]
        
        # If we have at least 4 characters to look back
        if len(result) >= 4:
            last_four = result[-4:]
            # Check if last four characters match any pattern
            if last_four == 'ABCD':
                result += 'A'
            elif last_four == 'BCDE':
                result += 'B'
            elif last_four == 'CDEA':
                result += 'C'
            elif last_four == 'DEAB':
                result += 'D'
            elif last_four == 'EABC':
                result += 'E'
        i += 1
    
    print(result)

s = "EABCBDEABADEAB"
apply_T(s)