def apply_t_operation(s):
    # Define the substrings and their corresponding insertions
    operations = {
        "ABCD": "A",
        "BCDE": "B",
        "CDEA": "C",
        "DEAB": "D",
        "EABC": "E"
    }
    
    # List to store the positions and characters to insert
    insertions = []
    
    # Iterate over the string to find all substrings
    i = 0
    while i <= len(s) - 4:
        substring = s[i:i+4]
        if substring in operations:
            # Record the position and the character to insert
            insertions.append((i + 4, operations[substring]))
        i += 1
    
    # Apply the insertions from left to right
    result = []
    last_pos = 0
    for pos, char in insertions:
        # Append the part of the string before the insertion
        result.append(s[last_pos:pos])
        # Append the character to insert
        result.append(char)
        # Update the last position
        last_pos = pos
    
    # Append the remaining part of the string
    result.append(s[last_pos:])
    
    # Join the result list into a final string
    return ''.join(result)

# Given string
s = "CADEABACDEADCDEAED"
# Apply the T operation
result = apply_t_operation(s)
print(result)