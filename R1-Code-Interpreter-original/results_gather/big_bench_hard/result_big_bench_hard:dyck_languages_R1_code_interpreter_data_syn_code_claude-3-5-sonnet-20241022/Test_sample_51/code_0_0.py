def find_missing_brackets(sequence):
    # Initialize stack for tracking brackets
    stack = []
    # Map of closing to opening brackets
    pairs = {')': '(', '}': '{', '>': '<'}
    # Map of opening to closing brackets
    opening_to_closing = {'(': ')', '{': '}', '<': '>'}
    
    # Process the given sequence
    for char in sequence:
        if char in '({<':
            stack.append(char)
        elif char in ')}>' and stack:
            if stack[-1] == pairs[char]:
                stack.pop()
            else:
                stack.append(char)
    
    # Generate the missing closing brackets
    missing = ''
    # Reverse the stack to get proper closing order
    stack.reverse()
    
    # Convert opening brackets to their closing counterparts
    for bracket in stack:
        if bracket in opening_to_closing:
            missing += opening_to_closing[bracket]
    
    print(f"Missing brackets: {missing}")

# Test the sequence
sequence = "(<{(){}}<>)>"
find_missing_brackets(sequence)