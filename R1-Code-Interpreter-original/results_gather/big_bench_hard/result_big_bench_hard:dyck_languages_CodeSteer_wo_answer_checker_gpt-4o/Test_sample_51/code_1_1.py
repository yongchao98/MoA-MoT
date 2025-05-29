def complete_sequence(sequence):
    # Define matching pairs
    matching_pairs = {'(': ')', '{': '}', '<': '>'}
    
    # Stack to keep track of opening symbols
    stack = []
    
    # Traverse the sequence
    for char in sequence:
        if char in matching_pairs:  # If it's an opening symbol
            stack.append(char)
        elif char in matching_pairs.values():  # If it's a closing symbol
            if stack and matching_pairs[stack[-1]] == char:
                stack.pop()
    
    # Complete the sequence by closing remaining symbols in the stack
    completion = ''.join(matching_pairs[char] for char in reversed(stack))
    
    # Output the final answer
    final_answer = sequence + completion
    print(f"<<<{final_answer}>>>")

# Given input
input_sequence = "( < { ( ) { } } ( < > ) >"
complete_sequence(input_sequence)