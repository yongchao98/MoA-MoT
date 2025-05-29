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
    offset = 0  # Keep track of how many characters we've inserted
    
    # Check each position in the original string
    for i in range(len(s)):
        # Get the substring of length 4 starting at position i
        if i + 4 <= len(s):
            substr = s[i:i+4]
            # Check if this substring matches any pattern
            if substr in patterns:
                # Insert the corresponding character at the correct position
                result.insert(i + 4 + offset, patterns[substr])
                offset += 1
    
    # Convert back to string and return
    return ''.join(result)

# Test with the given string
s = "EAABCDCABCDADEAB"
result = apply_T(s)
print(result)