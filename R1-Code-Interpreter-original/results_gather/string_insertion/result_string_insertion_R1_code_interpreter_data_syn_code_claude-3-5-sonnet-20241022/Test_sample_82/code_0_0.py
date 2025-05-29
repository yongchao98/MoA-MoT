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
    
    # Scan through the string
    for i in range(len(s)):
        # Get the substring starting at current position
        substr = s[i:i+4]
        # Check each rule
        for pattern, insert_char in rules:
            if substr == pattern:
                # Insert the character after the pattern
                # Account for previous insertions
                insert_pos = i + 4 + insertions
                result.insert(insert_pos, insert_char)
                insertions += 1
    
    return ''.join(result)

# Test with the given string
input_str = "ABCDBDEABAABCD"
result = apply_T(input_str)
print(result)