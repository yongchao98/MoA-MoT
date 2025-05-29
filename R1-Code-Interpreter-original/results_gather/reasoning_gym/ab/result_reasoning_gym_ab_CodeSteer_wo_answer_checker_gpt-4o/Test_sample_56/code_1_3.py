def compute_program(tokens):
    i = 0
    while i < len(tokens) - 1:
        if tokens[i] == 'A#' and tokens[i + 1] == '#A':
            # Remove both tokens
            tokens = tokens[:i] + tokens[i+2:]
            i = max(i - 1, 0)  # Step back to check the new pair
        elif tokens[i] == 'A#' and tokens[i + 1] == '#B':
            # Swap them to #B A#
            tokens[i], tokens[i + 1] = '#B', 'A#'
            i += 1
        elif tokens[i] == 'B#' and tokens[i + 1] == '#A':
            # Swap them to #A B#
            tokens[i], tokens[i + 1] = '#A', 'B#'
            i += 1
        elif tokens[i] == 'B#' and tokens[i + 1] == '#B':
            # Remove both tokens
            tokens = tokens[:i] + tokens[i+2:]
            i = max(i - 1, 0)  # Step back to check the new pair
        else:
            i += 1
        # Debug statement to trace the changes
        print(f"Current tokens: {tokens}")

    return tokens

# Initial program
program = ['#B', 'B#', '#A', 'B#', 'A#', 'A#', '#A', '#B', 'B#', 'B#']

# Compute the final state
final_state = compute_program(program)

# Output the final state
print(f"<<<{' '.join(final_state)}>>>")