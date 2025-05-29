def analyze_brackets(sequence):
    # Remove spaces
    sequence = sequence.replace(" ", "")
    
    # Stack to track brackets
    stack = []
    
    # Process sequence and show stack status
    print("Processing sequence:", sequence)
    print("\nStack progression:")
    
    for i, char in enumerate(sequence):
        if char in '[{<':
            stack.append(char)
            print(f"After {char}: {stack}")
        elif char in ']}>' and stack:
            opening = stack.pop()
            print(f"After {char}: {stack}")
    
    print("\nRemaining stack:", stack)
    if stack:
        print("Missing closures:", [']}>'['[{<'.index(x)] for x in stack])

analyze_brackets("[ < [ ] > ] [ { [ ] }")