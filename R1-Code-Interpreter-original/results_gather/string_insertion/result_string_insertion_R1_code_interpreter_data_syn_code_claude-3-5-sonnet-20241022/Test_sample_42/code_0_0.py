def apply_T(s):
    # Define the rules as dictionary: substring -> character to insert
    rules = {
        'ABCD': 'A',
        'BCDE': 'B',
        'CDEA': 'C',
        'DEAB': 'D',
        'EABC': 'E'
    }
    
    # Convert to list for easier manipulation
    result = list(s)
    # Keep track of insertions to adjust positions
    insertions = 0
    
    # Scan through the string
    for i in range(len(s)):
        # Check each rule
        for pattern, insert_char in rules.items():
            # Check if pattern starts at current position
            if i + len(pattern) <= len(s):
                substring = s[i:i + len(pattern)]
                if substring == pattern:
                    # Insert the character after the pattern
                    insert_pos = i + len(pattern) + insertions
                    result.insert(insert_pos, insert_char)
                    insertions += 1
    
    # Convert back to string
    return ''.join(result)

# Test with the given string
input_string = "EABCCEABCACDEA"
result = apply_T(input_string)
print(result)