def complete_sequence(sequence):
    # Define matching pairs
    matching_pairs = {']': '[', '}': '{', ')': '('}
    opening_chars = set(matching_pairs.values())
    
    # Initialize a stack
    stack = []
    
    # Iterate through the sequence
    for char in sequence:
        if char in opening_chars:
            stack.append(char)
        elif char in matching_pairs:
            if stack and stack[-1] == matching_pairs[char]:
                stack.pop()
            else:
                # Unbalanced sequence
                return "Unbalanced sequence"
    
    # Construct the final sequence
    closing_sequence = ''
    while stack:
        opening_char = stack.pop()
        for closing_char, opening_match in matching_pairs.items():
            if opening_match == opening_char:
                closing_sequence += closing_char
                break
    
    # Output the final answer
    return f"<<<{sequence}{closing_sequence}>>>"

# Initial sequence
initial_sequence = "[ [ [ [ < { ( ) } ( ) > { ( ) } ] { ( ) }"

# Complete the sequence
final_answer = complete_sequence(initial_sequence)

# Print the final answer
print(final_answer)