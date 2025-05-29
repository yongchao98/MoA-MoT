def balance_parentheses(sequence):
    open_count = 0
    for char in sequence:
        if char == '(':
            open_count += 1
        elif char == ')':
            if open_count > 0:
                open_count -= 1
    return ')' * open_count

sequence = "(((()"
closing_needed = balance_parentheses(sequence)
print(sequence + closing_needed)