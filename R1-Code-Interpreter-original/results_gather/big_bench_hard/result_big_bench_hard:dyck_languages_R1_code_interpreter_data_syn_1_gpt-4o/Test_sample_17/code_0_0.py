def complete_parentheses(sequence):
    balance = 0
    for char in sequence:
        if char == '(':
            balance += 1
        elif char == ')':
            balance -= 1
    
    # Add the necessary number of closing parentheses
    completed_sequence = sequence + ')' * balance
    return completed_sequence

sequence = "< ( ( )"
completed_sequence = complete_parentheses(sequence)
print(completed_sequence)