def is_balanced(sequence):
    stack = []
    brackets = {'(': ')', '[': ']', '<': '>'}
    
    for char in sequence:
        if char in '([<':
            stack.append(char)
        elif char in ')]>':
            if not stack:
                return False
            if char != brackets[stack.pop()]:
                return False
    
    return len(stack) == 0

# Test original sequence and completed sequence
original = "<<[()]>"
completed = "<<[()]>>"

print(f"Original sequence '{original}' is balanced: {is_balanced(original)}")
print(f"Completed sequence '{completed}' is balanced: {is_balanced(completed)}")