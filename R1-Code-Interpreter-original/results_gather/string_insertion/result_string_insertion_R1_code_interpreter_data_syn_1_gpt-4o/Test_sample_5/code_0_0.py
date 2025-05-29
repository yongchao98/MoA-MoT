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
    i = 0
    while i <= len(s) - 4:
        for sub, char in operations:
            if s[i:i+4] == sub:
                insertions.append((i+4, char))
        i += 1
    
    # Sort insertions by position (though they should already be in order)
    insertions.sort()
    
    # Create a new string with the insertions
    result = []
    last_pos = 0
    for pos, char in insertions:
        result.append(s[last_pos:pos])
        result.append(char)
        last_pos = pos
    result.append(s[last_pos:])
    
    return ''.join(result)

# Given string
s = "DEABBEABCCABCDA"

# Apply the T operation
result = apply_t_operation(s)

# Output the result
print(result)