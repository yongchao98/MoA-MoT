def complete_sequence(input_sequence):
    open_count = 0
    for char in input_sequence:
        if char == '(':
            open_count += 1
        elif char == ')':
            open_count -= 1

    # Add the necessary number of closing parentheses
    completed_sequence = input_sequence + ')' * open_count
    return f"<<<{completed_sequence}>>>"

input_sequence = "< ( ( )"
print(complete_sequence(input_sequence))