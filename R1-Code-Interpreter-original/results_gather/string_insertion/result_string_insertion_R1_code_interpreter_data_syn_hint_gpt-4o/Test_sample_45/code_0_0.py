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
            insertions.append((i+4, operations[substring]))
        i += 1
    
    # Apply insertions from left to right
    result = []
    last_pos = 0
    for pos, char in insertions:
        result.append(s[last_pos:pos])
        result.append(char)
        last_pos = pos
    result.append(s[last_pos:])
    
    return ''.join(result)

# Given string
s = "DDBCDECABCDCCDEA"
# Apply the T operation
result = apply_t_operation(s)
print(result)