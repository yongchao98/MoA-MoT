def apply_T(s):
    # Define the rules as tuples (pattern, insert_char)
    rules = [
        ('ABCD', 'A'),
        ('BCDE', 'B'),
        ('CDEA', 'C'),
        ('DEAB', 'D'),
        ('EABC', 'E')
    ]
    
    # Convert to list for easier manipulation
    result = list(s)
    # Keep track of insertions to adjust positions
    insertions = 0
    
    # Check each position in the original string
    for i in range(len(s)):
        # Get the substring starting at current position
        # We need 4 characters for checking
        if i + 4 <= len(s):
            substr = s[i:i+4]
            # Check against each rule
            for pattern, insert_char in rules:
                if substr == pattern:
                    # Insert after the pattern, accounting for previous insertions
                    result.insert(i + 4 + insertions, insert_char)
                    insertions += 1
    
    # Convert back to string
    return ''.join(result)

# Test with the given string
input_str = "ABCDBDEABEEABC"
result = apply_T(input_str)
print(result)