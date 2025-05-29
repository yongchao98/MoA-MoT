def apply_t_operation(s):
    # Define the substrings and their corresponding insertions
    operations = [
        ("ABCD", "A"),
        ("BCDE", "B"),
        ("CDEA", "C"),
        ("DEAB", "D"),
        ("EABC", "E")
    ]
    
    # List to store the positions and characters to insert
    insertions = []
    
    # Iterate over the string to find all substrings
    for i in range(len(s) - 3):
        for sub, char in operations:
            if s[i:i+4] == sub:
                insertions.append((i + 4, char))
    
    # Sort insertions by position (though they should already be in order)
    insertions.sort()
    
    # Create a new list to build the result string
    result = []
    last_pos = 0
    
    # Insert characters into the result string
    for pos, char in insertions:
        result.append(s[last_pos:pos])
        result.append(char)
        last_pos = pos
    
    # Append the remaining part of the original string
    result.append(s[last_pos:])
    
    # Join the result list into a final string
    return ''.join(result)

# Given string
s = "EABCBABCDADEABDB"

# Apply the T operation
result = apply_t_operation(s)

# Output the result
print(result)