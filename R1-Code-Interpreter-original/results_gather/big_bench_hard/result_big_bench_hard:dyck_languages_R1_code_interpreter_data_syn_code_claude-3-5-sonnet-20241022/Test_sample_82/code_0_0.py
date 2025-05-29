def find_missing_brackets(sequence):
    # Define matching pairs
    pairs = {
        '<': '>', 
        '{': '}', 
        '(': ')'
    }
    
    # Remove spaces and convert to list
    sequence = ''.join(sequence.split())
    
    # Stack to track opening brackets
    stack = []
    
    # Process existing sequence
    for char in sequence:
        if char in '<{(':  # Opening brackets
            stack.append(char)
        elif char in '>})':  # Closing brackets
            if not stack:
                continue
            if (char == '>' and stack[-1] == '<') or \
               (char == '}' and stack[-1] == '{') or \
               (char == ')' and stack[-1] == '('):
                stack.pop()
    
    # Generate missing closing brackets
    missing = ''
    while stack:
        missing = pairs[stack.pop()] + missing
    
    print(f"Original sequence: {sequence}")
    print(f"Missing brackets: {missing}")
    print(f"Complete sequence: {sequence + missing}")

# Test the sequence
find_missing_brackets("< { ( { } ) < > }")