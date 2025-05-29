def complete_sequence(sequence):
    # Define matching pairs
    matching_pairs = {'(': ')', '{': '}', '[': ']', '<': '>'}
    
    # Stack to keep track of opening brackets
    stack = []
    
    # Traverse the sequence
    for char in sequence:
        if char in matching_pairs:  # If it's an opening bracket
            stack.append(char)
        elif char in matching_pairs.values():  # If it's a closing bracket
            if stack and matching_pairs[stack[-1]] == char:
                stack.pop()
            else:
                return "Invalid sequence"
    
    # If stack is not empty, add the necessary closing brackets
    completion = ''.join(matching_pairs[char] for char in reversed(stack))
    
    # Output the final answer
    print(f"<<<{sequence}{completion}>>>")

# Input sequence
input_sequence = "( ( { < { < { < > } > } > } ) ( ) ( { [ ( < ( < < { ( { < > } ) } > > ) > ) { ( { < ( { [ ] } ( ) ) > } < { < { { ( ) } } > { < [ { ( ) } ] > } } > ) } [ [ ( ) ] ] ]"

# Complete the sequence
complete_sequence(input_sequence)