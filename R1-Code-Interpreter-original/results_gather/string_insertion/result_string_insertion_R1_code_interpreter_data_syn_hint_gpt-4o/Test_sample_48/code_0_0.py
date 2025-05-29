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
    
    # Find all occurrences of each substring and prepare insertions
    for substring, char_to_insert in operations:
        start = 0
        while True:
            index = s.find(substring, start)
            if index == -1:
                break
            # Record the position for insertion (after the substring)
            insertions.append((index + len(substring), char_to_insert))
            start = index + 1
    
    # Sort insertions by position to ensure left-to-right order
    insertions.sort()
    
    # Apply insertions to the string
    result = []
    last_index = 0
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
s = "EADEABDABCDCBCDE"

# Apply the T operation
resulting_string = apply_t_operation(s)

# Output the resulting string
print(resulting_string)