def compute_program(tokens):
    changed = True
    while changed:
        changed = False
        i = 0
        while i < len(tokens) - 1:
            if tokens[i] == 'A#' and tokens[i + 1] == '#A':
                # Apply rule: A# #A -> nothing
                tokens.pop(i)
                tokens.pop(i)
                changed = True
            elif tokens[i] == 'A#' and tokens[i + 1] == '#B':
                # Apply rule: A# #B -> #B A#
                tokens[i], tokens[i + 1] = '#B', 'A#'
                changed = True
                i += 1
            elif tokens[i] == 'B#' and tokens[i + 1] == '#A':
                # Apply rule: B# #A -> #A B#
                tokens[i], tokens[i + 1] = '#A', 'B#'
                changed = True
                i += 1
            elif tokens[i] == 'B#' and tokens[i + 1] == '#B':
                # Apply rule: B# #B -> nothing
                tokens.pop(i)
                tokens.pop(i)
                changed = True
            else:
                i += 1
            # If a change was made, reset the index to start from the beginning
            if changed:
                i = 0
    return tokens

# Initial program
program = ['#B', '#A', '#B', '#B', 'A#', '#A', 'B#', '#A', 'B#', '#A']

# Compute the final state
final_state = compute_program(program)

# Output the final state
print(f"<<<{' '.join(final_state)}>>>")