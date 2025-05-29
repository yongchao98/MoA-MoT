def analyze_nesting(sequence):
    stack = []
    nesting_level = 0
    pairs = {')': '(', '}': '{', '>': '<'}
    opening_to_closing = {'(': ')', '{': '}', '<': '>'}
    
    print("Analyzing sequence:", sequence)
    print("\nStep by step analysis:")
    
    for i, char in enumerate(sequence):
        print(f"\nPosition {i}, Character: {char}")
        if char in '({<':
            nesting_level += 1
            stack.append((char, nesting_level))
            print(f"Opening at level {nesting_level}")
            print(f"Stack: {[f'{c}(L{l})' for c, l in stack]}")
        elif char in ')}>' and stack:
            expected_open = pairs[char]
            if stack and stack[-1][0] == expected_open:
                open_char, level = stack.pop()
                nesting_level = level - 1
                print(f"Matched with {open_char} from level {level}")
            else:
                print(f"Mismatch: Expected {pairs[char]}")
            print(f"Stack: {[f'{c}(L{l})' for c, l in stack]}")
    
    print("\nRemaining brackets:")
    for char, level in stack:
        print(f"{char} at level {level} needs {opening_to_closing[char]}")

sequence = "(<{(){}}<>)>"
analyze_nesting(sequence)