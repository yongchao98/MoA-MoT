def find_closing_sequence(input_str):
    # Define matching brackets
    brackets = {
        '(': ')',
        '[': ']',
        '{': '}',
        '<': '>'
    }
    
    # Initialize stack to keep track of open brackets
    stack = []
    
    # Process input string
    for char in input_str:
        if char in brackets.keys():  # Opening bracket
            stack.append(char)
        elif char in brackets.values():  # Closing bracket
            if stack and brackets[stack[-1]] == char:
                stack.pop()
            else:
                return "Invalid sequence"
    
    # Generate closing sequence
    closing_sequence = ''
    while stack:
        closing_sequence = brackets[stack.pop()] + closing_sequence
    
    # Print original sequence and completion
    print(f"Original sequence: {input_str}")
    print(f"Required closing brackets: {closing_sequence}")
    print(f"Complete sequence: {input_str + closing_sequence}")

# Test with given input
input_str = "[ ( < { [ ( ) ] } >"
find_closing_sequence(input_str)