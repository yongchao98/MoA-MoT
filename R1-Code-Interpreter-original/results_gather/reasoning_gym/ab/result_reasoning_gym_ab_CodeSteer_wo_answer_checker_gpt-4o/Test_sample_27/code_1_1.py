def compute_program(tokens):
    stack = []
    for token in tokens:
        if stack:
            top = stack[-1]
            if top == 'A#' and token == '#A':
                stack.pop()  # A# #A becomes nothing
            elif top == 'B#' and token == '#B':
                stack.pop()  # B# #B becomes nothing
            elif top == 'A#' and token == '#B':
                stack[-1] = '#B'  # A# #B becomes #B A#
                stack.append('A#')
            elif top == 'B#' and token == '#A':
                stack[-1] = '#A'  # B# #A becomes #A B#
                stack.append('B#')
            else:
                stack.append(token)
        else:
            stack.append(token)
    return stack

# Initial program
program = ['B#', 'A#', '#B', 'A#', 'B#', '#A', 'A#', 'A#', '#B', '#B']

# Compute the final state
final_state = compute_program(program)

# Output the final state
print(f"<<<{' '.join(final_state)}>>>")