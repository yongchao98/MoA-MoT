def compute_program(tokens):
    stack = []
    for token in tokens:
        if stack:
            top = stack[-1]
            if (top == 'A#' and token == '#A') or (top == 'B#' and token == '#B'):
                # Remove both tokens
                stack.pop()
            elif top == 'A#' and token == '#B':
                # Swap them
                stack.pop()
                stack.append('#B')
                stack.append('A#')
            elif top == 'B#' and token == '#A':
                # Swap them
                stack.pop()
                stack.append('#A')
                stack.append('B#')
            else:
                # No rule applies, push the current token
                stack.append(token)
        else:
            # Stack is empty, push the current token
            stack.append(token)
    return stack

# Initial program
program = ['#A', '#B', '#B', '#A', 'A#', '#A', '#A', 'A#', 'A#', '#B']

# Compute the final state
final_state = compute_program(program)

# Output the final state
print(f"<<<{' '.join(final_state)}>>>")