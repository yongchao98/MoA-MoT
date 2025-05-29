def compute_program(tokens):
    while True:
        i = 0
        changes_made = False
        while i < len(tokens) - 1:
            if tokens[i] == 'A#' and tokens[i + 1] == '#A':
                # Apply rule: A# #A becomes nothing
                tokens.pop(i)
                tokens.pop(i)
                changes_made = True
            elif tokens[i] == 'A#' and tokens[i + 1] == '#B':
                # Apply rule: A# #B becomes #B A#
                tokens[i], tokens[i + 1] = '#B', 'A#'
                changes_made = True
                i += 1
            elif tokens[i] == 'B#' and tokens[i + 1] == '#A':
                # Apply rule: B# #A becomes #A B#
                tokens[i], tokens[i + 1] = '#A', 'B#'
                changes_made = True
                i += 1
            elif tokens[i] == 'B#' and tokens[i + 1] == '#B':
                # Apply rule: B# #B becomes nothing
                tokens.pop(i)
                tokens.pop(i)
                changes_made = True
            else:
                i += 1
        
        if not changes_made:
            break

    return tokens

# Initial program
program = ['#B', 'A#', '#B', 'A#', 'B#', '#A', 'B#', 'B#', 'A#', '#B']

# Compute the final state
final_state = compute_program(program)

# Output the final state
print(f"<<<{' '.join(final_state)}>>>")