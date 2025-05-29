def analyze_brackets(s):
    # Define matching pairs
    pairs = {
        '(': ')',
        '[': ']',
        '{': '}',
        '<': '>'
    }
    
    # Initialize stack and counter
    stack = []
    count = {')': 0, ']': 0, '}': 0, '>': 0}
    opening_count = {'(': 0, '[': 0, '{': 0, '<': 0}
    
    # Process the input string
    for char in s:
        if char in '([{<':
            stack.append(char)
            opening_count[char] += 1
        elif char in ')]}>' and stack:
            if pairs[stack[-1]] == char:
                count[char] += 1
                stack.pop()
            else:
                print("Invalid nesting")
                return
    
    # Print detailed analysis
    print("Opening brackets:", opening_count)
    print("Closing brackets:", count)
    print("Remaining stack:", ''.join(stack))
    
    # Generate missing brackets
    missing = ''
    while stack:
        missing = pairs[stack.pop()] + missing
    print("Missing brackets needed:", missing)

# Test with input
input_seq = "< [ ] > ( ( ( { { { { ( ( { { } } ) ) } ( { } ) } } { { { [ { ( ) } ( ) ( ( { [ < < > > ] < < > > } ) ) ] ( { } { } ) } } } } )"
analyze_brackets(input_seq)