def analyze_brackets(sequence):
    stack = []
    pairs = {')': '(', '}': '{', '>': '<'}
    opening = '({<'
    closing = ')}>'
    
    print("Processing sequence:", sequence)
    print("\nStep by step analysis:")
    
    for i, char in enumerate(sequence):
        print(f"\nPosition {i}, Character: {char}")
        if char in opening:
            stack.append(char)
            print(f"Added opening bracket: {char}")
            print(f"Current stack: {stack}")
        elif char in closing:
            if stack and stack[-1] == pairs[char]:
                removed = stack.pop()
                print(f"Matched and removed: {removed}")
                print(f"Current stack: {stack}")
            else:
                print(f"Mismatch or empty stack")
                print(f"Current stack: {stack}")
    
    print("\nFinal unmatched brackets:", stack)
    
    # Generate missing closing brackets
    missing = []
    for bracket in reversed(stack):
        if bracket == '(':
            missing.append(')')
        elif bracket == '{':
            missing.append('}')
        elif bracket == '<':
            missing.append('>')
    
    print("Missing closing brackets (in order):", missing)

sequence = "(<{(){}}<>)>"
analyze_brackets(sequence)