def compute_program(tokens):
    stack = []
    i = 0
    while i < len(tokens):
        if stack:
            top = stack[-1]
            current = tokens[i]
            if top == 'A#' and current == '#A':
                stack.pop()  # A# #A becomes nothing
            elif top == 'A#' and current == '#B':
                stack[-1] = '#B'  # A# #B becomes #B A#
                stack.append('A#')
            elif top == 'B#' and current == '#A':
                stack[-1] = '#A'  # B# #A becomes #A B#
                stack.append('B#')
            elif top == 'B#' and current == '#B':
                stack.pop()  # B# #B becomes nothing
            else:
                stack.append(current)
        else:
            stack.append(tokens[i])
        i += 1
    return stack

initial_program = ['B#', 'B#', 'A#', '#A', '#A', '#A', 'A#', '#A', 'A#', '#A']
final_program = compute_program(initial_program)
print(' '.join(final_program))