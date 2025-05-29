def compute_program(tokens):
    while True:
        i = 0
        changed = False
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
            # Debugging statement to track changes
            print(f"Current state: {' '.join(tokens)}")
        if not changed:
            break
    return tokens

# Initial program
program = ['#A', 'B#', '#B', 'A#', 'A#', '#B', 'B#', 'B#', '#A', '#A']

# Compute the final state
final_state = compute_program(program)

# Output the final state
print(f"<<<{' '.join(final_state)}>>>")