def compute_program(tokens):
    i = 0
    while i < len(tokens) - 1:
        if tokens[i] == 'A#' and tokens[i + 1] == '#A':
            # Apply A# #A -> nothing
            tokens.pop(i)
            tokens.pop(i)
            i = max(i - 1, 0)  # Move back one step to recheck
        elif tokens[i] == 'A#' and tokens[i + 1] == '#B':
            # Apply A# #B -> #B A#
            tokens[i], tokens[i + 1] = tokens[i + 1], tokens[i]
            i += 1
        elif tokens[i] == 'B#' and tokens[i + 1] == '#A':
            # Apply B# #A -> #A B#
            tokens[i], tokens[i + 1] = tokens[i + 1], tokens[i]
            i += 1
        elif tokens[i] == 'B#' and tokens[i + 1] == '#B':
            # Apply B# #B -> nothing
            tokens.pop(i)
            tokens.pop(i)
            i = max(i - 1, 0)  # Move back one step to recheck
        else:
            i += 1
    return tokens

# Initial program
program = ['A#', '#A', 'A#', 'B#', '#A', '#B', '#A', '#B', '#B', 'A#']
final_state = compute_program(program)
print(final_state)