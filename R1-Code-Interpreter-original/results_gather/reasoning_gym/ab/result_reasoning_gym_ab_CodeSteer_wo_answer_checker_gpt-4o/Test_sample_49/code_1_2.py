def compute_program(tokens):
    changed = True
    while changed:
        changed = False
        i = 0
        while i < len(tokens) - 1:
            if tokens[i] == 'A#' and tokens[i + 1] == '#A':
                # Remove both tokens
                tokens.pop(i)
                tokens.pop(i)
                changed = True
                i = max(i - 1, 0)  # Step back to check new neighbors
            elif tokens[i] == 'A#' and tokens[i + 1] == '#B':
                # Replace with #B A#
                tokens[i], tokens[i + 1] = '#B', 'A#'
                changed = True
                i += 1
            elif tokens[i] == 'B#' and tokens[i + 1] == '#A':
                # Replace with #A B#
                tokens[i], tokens[i + 1] = '#A', 'B#'
                changed = True
                i += 1
            elif tokens[i] == 'B#' and tokens[i + 1] == '#B':
                # Remove both tokens
                tokens.pop(i)
                tokens.pop(i)
                changed = True
                i = max(i - 1, 0)  # Step back to check new neighbors
            else:
                i += 1
    return tokens

# Initial program
program = ['#B', '#A', '#B', '#B', 'A#', 'B#', 'B#', '#B', 'A#', '#A']

# Compute the final state
final_state = compute_program(program)

# Output the final state
print(f"<<<{' '.join(final_state)}>>>")