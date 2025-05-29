def compute_program(tokens):
    i = 0
    while i < len(tokens) - 1:
        print(f"Current tokens: {tokens}")  # Debug statement
        if tokens[i] == 'A#' and tokens[i + 1] == '#A':
            # Remove both tokens
            tokens.pop(i)
            tokens.pop(i)
            i = max(i - 1, 0)  # Move back one step to recheck
        elif tokens[i] == 'B#' and tokens[i + 1] == '#B':
            # Remove both tokens
            tokens.pop(i)
            tokens.pop(i)
            i = max(i - 1, 0)  # Move back one step to recheck
        elif tokens[i] == 'A#' and tokens[i + 1] == '#B':
            # Swap them
            tokens[i], tokens[i + 1] = '#B', 'A#'
            i += 1
        elif tokens[i] == 'B#' and tokens[i + 1] == '#A':
            # Swap them
            tokens[i], tokens[i + 1] = '#A', 'B#'
            i += 1
        else:
            i += 1
    print(f"Final tokens: {tokens}")  # Debug statement
    return tokens

# Initial program
program = ['#A', '#B', '#B', '#A', 'A#', '#A', '#A', 'A#', 'A#', '#B']

# Compute the final state
final_state = compute_program(program)

# Output the final state
print(f"<<<{' '.join(final_state)}>>>")