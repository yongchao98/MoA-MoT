def find_missing_brackets(sequence):
    # Define matching pairs
    pairs = {
        '[': ']',
        '{': '}',
        '<': '>',
        '(': ')'
    }
    
    # Convert input to list for easier manipulation
    brackets = list(sequence.replace(" ", ""))
    
    # Stack to keep track of opening brackets
    stack = []
    
    # Process existing sequence
    for bracket in brackets:
        if bracket in pairs.keys():  # Opening bracket
            stack.append(bracket)
        elif bracket in pairs.values():  # Closing bracket
            if not stack:
                continue
            if bracket == pairs[stack[-1]]:
                stack.pop()
            
    # Generate missing closing brackets
    missing = ''
    while stack:
        missing = pairs[stack.pop()] + missing
        
    # Print original sequence and completed sequence
    print(f"Original sequence: {sequence}")
    print(f"Missing brackets: {missing}")
    print(f"Complete sequence: {sequence + missing}")

# Test the function
find_missing_brackets("[ { < ( ) > }")