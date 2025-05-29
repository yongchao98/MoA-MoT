def apply_t_operation(s):
    # Define the operations as a list of tuples (substring, character to insert)
    operations = [
        ("ABCD", "A"),
        ("BCDE", "B"),
        ("CDEA", "C"),
        ("DEAB", "D"),
        ("EABC", "E")
    ]
    
    # Convert the string to a list for easier manipulation
    s_list = list(s)
    i = 0
    
    # Iterate over the string
    while i < len(s_list) - 3:
        # Check each operation
        for sub, char in operations:
            # If the substring is found, insert the character
            if s_list[i:i+4] == list(sub):
                s_list.insert(i+4, char)
                # Move the index past the inserted character
                i += 4
                break
        i += 1
    
    # Convert the list back to a string
    return ''.join(s_list)

# Given string
s = "CCDEAACDEACEABCBD"
# Apply the T operation
result = apply_t_operation(s)
print(result)