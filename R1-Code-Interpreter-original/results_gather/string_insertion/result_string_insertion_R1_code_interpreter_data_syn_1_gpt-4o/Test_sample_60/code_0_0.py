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
    
    # Find all occurrences of each substring and store the insertion positions
    for substring, char_to_insert in operations:
        start = 0
        while True:
            index = s.find(substring, start)
            if index == -1:
                break
            # Append the position where the character should be inserted
            insertions.append((index + len(substring), char_to_insert))
            start = index + 1
    
    # Sort insertions by position to ensure left-to-right order
    insertions.sort()
    
    # Create a new list to build the result string
    result = []
    last_index = 0
    
    # Insert characters into the result string
    for index, char in insertions:
        # Append the part of the string before the insertion point
        result.append(s[last_index:index])
        # Append the character to insert
        result.append(char)
        # Update the last index
        last_index = index
    
    # Append the remaining part of the string
    result.append(s[last_index:])
    
    # Join the result list into a final string
    return ''.join(result)

# Given string
s = "BCDEDDEABBEABCAE"

# Apply the T operation
result = apply_t_operation(s)

# Print the result
print(result)