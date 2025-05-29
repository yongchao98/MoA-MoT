def compute_program(tokens):
    i = 0
    while i < len(tokens) - 1:
        if tokens[i] == 'A#' and tokens[i + 1] == '#A':
            tokens.pop(i)
            tokens.pop(i)
            i = 0  # Restart from the beginning
        elif tokens[i] == 'A#' and tokens[i + 1] == '#B':
            tokens[i], tokens[i + 1] = '#B', 'A#'
            i = 0  # Restart from the beginning
        elif tokens[i] == 'B#' and tokens[i + 1] == '#A':
            tokens[i], tokens[i + 1] = '#A', 'B#'
            i = 0  # Restart from the beginning
        elif tokens[i] == 'B#' and tokens[i + 1] == '#B':
            tokens.pop(i)
            tokens.pop(i)
            i = 0  # Restart from the beginning
        else:
            i += 1
    return tokens

# Initial program
program = ['#B', '#B', '#B', 'A#', '#B', '#A', 'A#', '#A', 'A#', '#A']

# Compute the final state
final_state = compute_program(program)

# Print the final state
print(" ".join(final_state))