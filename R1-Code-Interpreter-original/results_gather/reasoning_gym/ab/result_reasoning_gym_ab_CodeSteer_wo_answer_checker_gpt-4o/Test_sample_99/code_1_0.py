def compute_program(tokens):
    while True:
        new_tokens = []
        i = 0
        changed = False
        while i < len(tokens):
            if i < len(tokens) - 1:
                if tokens[i] == 'A#' and tokens[i + 1] == '#A':
                    # Apply A# #A -> nothing
                    i += 2
                    changed = True
                elif tokens[i] == 'A#' and tokens[i + 1] == '#B':
                    # Apply A# #B -> #B A#
                    new_tokens.extend(['#B', 'A#'])
                    i += 2
                    changed = True
                elif tokens[i] == 'B#' and tokens[i + 1] == '#A':
                    # Apply B# #A -> #A B#
                    new_tokens.extend(['#A', 'B#'])
                    i += 2
                    changed = True
                elif tokens[i] == 'B#' and tokens[i + 1] == '#B':
                    # Apply B# #B -> nothing
                    i += 2
                    changed = True
                else:
                    new_tokens.append(tokens[i])
                    i += 1
            else:
                new_tokens.append(tokens[i])
                i += 1

        if not changed:
            break
        tokens = new_tokens

    return tokens

# Initial program
program = ['#A', '#A', 'A#', '#B', '#A', 'A#', '#B', 'A#', '#A', 'A#']

# Compute the final state
final_state = compute_program(program)

# Output the final state
print(f"<<<{' '.join(final_state)}>>>")