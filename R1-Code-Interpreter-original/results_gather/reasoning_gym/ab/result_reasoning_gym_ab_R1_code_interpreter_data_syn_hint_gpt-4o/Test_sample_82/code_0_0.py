def compute_program(tokens):
    i = 0
    while i < len(tokens) - 1:
        if tokens[i] == 'A#' and tokens[i + 1] == '#A':
            # A# #A becomes nothing
            tokens.pop(i)
            tokens.pop(i)
            i = 0  # Restart from the beginning
        elif tokens[i] == 'A#' and tokens[i + 1] == '#B':
            # A# #B becomes #B A#
            tokens[i], tokens[i + 1] = '#B', 'A#'
            i = 0  # Restart from the beginning
        elif tokens[i] == 'B#' and tokens[i + 1] == '#A':
            # B# #A becomes #A B#
            tokens[i], tokens[i + 1] = '#A', 'B#'
            i = 0  # Restart from the beginning
        elif tokens[i] == 'B#' and tokens[i + 1] == '#B':
            # B# #B becomes nothing
            tokens.pop(i)
            tokens.pop(i)
            i = 0  # Restart from the beginning
        else:
            i += 1  # Move to the next pair
    return tokens

# Initial program
program = ['B#', 'A#', '#A', 'B#', 'A#', 'B#', '#A', 'A#', '#A', '#A']

# Compute the final state of the program
final_state = compute_program(program)

# Print the final state
print(final_state)