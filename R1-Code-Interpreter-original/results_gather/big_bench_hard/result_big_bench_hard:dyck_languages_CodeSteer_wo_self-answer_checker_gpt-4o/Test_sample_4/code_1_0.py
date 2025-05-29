def is_valid_sequence(sequence):
    # Define matching pairs
    matching_pairs = {'>': '<', ')': '(', '}': '{', ']': '['}
    # Initialize an empty stack
    stack = []
    
    for char in sequence:
        if char in '<({[':
            # Push opening brackets onto the stack
            stack.append(char)
            print(f"Pushed {char}, stack: {stack}")
        elif char in '>)}]':
            # Check if the stack is not empty and the top of the stack matches
            if stack and stack[-1] == matching_pairs[char]:
                stack.pop()
                print(f"Popped {char}, stack: {stack}")
            else:
                print(f"Mismatch or empty stack for {char}, stack: {stack}")
                return False
        else:
            # Ignore non-bracket characters
            continue
    
    # If the stack is empty, all brackets are properly closed
    if len(stack) == 0:
        return True
    else:
        print(f"Unmatched opening brackets remain, stack: {stack}")
        return False

sequence = "< ( ( ( < > ) ) ( { { } [ { } ] [ ] < ( ) > } ) )"
if is_valid_sequence(sequence):
    print(f"<<<{sequence}>>>")
else:
    print("The sequence is not properly closed.")