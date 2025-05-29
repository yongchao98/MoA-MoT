def apply_t_operation(s):
    # Define the substrings and their corresponding insertions
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
        for sub, insert in operations:
            # If the substring is found, perform the insertion
            if s_list[i:i+4] == list(sub):
                s_list.insert(i+4, insert)
                i += 4  # Move past the inserted character
                break
        else:
            i += 1  # Move to the next character if no insertion was made
    
    # Convert the list back to a string
    return ''.join(s_list)

# Given string
initial_string = "ADEABEDEABBEABCBE"
# Apply the T operation
result_string = apply_t_operation(initial_string)
print(result_string)