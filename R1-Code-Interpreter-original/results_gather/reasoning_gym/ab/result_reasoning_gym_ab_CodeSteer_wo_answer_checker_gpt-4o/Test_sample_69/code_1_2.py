def compute_program(tokens):
    changed = True
    iteration = 0
    while changed:
        changed = False
        i = 0
        print(f"Iteration {iteration}: {' '.join(tokens)}")  # Debugging output
        while i < len(tokens) - 1:
            if tokens[i] == 'A#' and tokens[i + 1] == '#A':
                # A# #A becomes nothing
                tokens.pop(i)
                tokens.pop(i)
                changed = True
            elif tokens[i] == 'A#' and tokens[i + 1] == '#B':
                # A# #B becomes #B A#
                tokens[i], tokens[i + 1] = '#B', 'A#'
                changed = True
                i += 1
            elif tokens[i] == 'B#' and tokens[i + 1] == '#A':
                # B# #A becomes #A B#
                tokens[i], tokens[i + 1] = '#A', 'B#'
                changed = True
                i += 1
            elif tokens[i] == 'B#' and tokens[i + 1] == '#B':
                # B# #B becomes nothing
                tokens.pop(i)
                tokens.pop(i)
                changed = True
            else:
                i += 1
        iteration += 1
    print(f"Final state: {' '.join(tokens)}")  # Debugging output
    return tokens

# Initial program
program = ['B#', '#A', '#B', 'A#', 'B#', 'B#', '#A', '#B', 'B#', '#A']

# Compute the final state
final_state = compute_program(program)

# Output the final state
print(f"<<<{' '.join(final_state)}>>>")