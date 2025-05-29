def complete_sequence(sequence):
    # Define matching pairs
    matching_bracket = {'<': '>', '[': ']', '{': '}', '(': ')'}
    
    # Stack to keep track of opening brackets
    stack = []
    
    # Traverse the sequence
    for char in sequence:
        if char in matching_bracket:  # If it's an opening bracket
            stack.append(char)
        elif char in matching_bracket.values():  # If it's a closing bracket
            if stack and matching_bracket[stack[-1]] == char:
                stack.pop()
    
    # Complete the sequence by adding the necessary closing brackets
    completion = ''.join(matching_bracket[char] for char in reversed(stack))
    
    # Output the final answer
    print(f"<<<{completion}>>>")

# Given input sequence
input_sequence = "[ < < > > < ( < < < [ [ [ [ ] ] ] ] > > > ) [ < < [ ] > > ] { [ < > ] < > [ ( ) ] } { [ ]"
complete_sequence(input_sequence)