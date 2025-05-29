def compute_program(tokens):
    stack = []
    i = 0
    while i < len(tokens):
        if stack:
            top = stack[-1]
            current = tokens[i]
            if (top == 'A#' and current == '#A') or (top == 'B#' and current == '#B'):
                # Remove both tokens
                stack.pop()
            elif top == 'A#' and current == '#B':
                # Swap them
                stack.pop()
                stack.append('#B')
                stack.append('A#')
            elif top == 'B#' and current == '#A':
                # Swap them
                stack.pop()
                stack.append('#A')
                stack.append('B#')
            else:
                # No rule applies, push the current token
                stack.append(current)
            i += 1
        else:
            # Stack is empty, push the current token
            stack.append(tokens[i])
            i += 1
    return stack

# Initial program
program = ['#A', '#B', '#B', '#A', 'A#', '#A', '#A', 'A#', 'A#', '#B']

# Compute the final state
final_state = compute_program(program)

# Output the final state
print(f"<<<{' '.join(final_state)}>>>")