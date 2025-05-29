def compute_program(tokens):
    i = 0
    while i < len(tokens) - 1:
        if tokens[i] == 'A#' and tokens[i + 1] == '#A':
            # A# #A becomes nothing
            tokens.pop(i)
            tokens.pop(i)
            i = max(i - 1, 0)  # Step back to check new neighbors
        elif tokens[i] == 'A#' and tokens[i + 1] == '#B':
            # A# #B becomes #B A#
            tokens[i], tokens[i + 1] = '#B', 'A#'
            i += 1
        elif tokens[i] == 'B#' and tokens[i + 1] == '#A':
            # B# #A becomes #A B#
            tokens[i], tokens[i + 1] = '#A', 'B#'
            i += 1
        elif tokens[i] == 'B#' and tokens[i + 1] == '#B':
            # B# #B becomes nothing
            tokens.pop(i)
            tokens.pop(i)
            i = max(i - 1, 0)  # Step back to check new neighbors
        else:
            i += 1
    return tokens

# Initial program
program = "#A A# A# #A #A #A B# #A B# B#".split()

# Compute the final state
final_state = compute_program(program)

# Output the final state
print(f"<<<{' '.join(final_state)}>>>")