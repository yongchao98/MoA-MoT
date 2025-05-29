def compute_program(tokens):
    while True:
        changed = False
        i = 0
        while i < len(tokens) - 1:
            if tokens[i] == 'A#' and tokens[i + 1] == '#A':
                tokens.pop(i)
                tokens.pop(i)
                changed = True
            elif tokens[i] == 'A#' and tokens[i + 1] == '#B':
                tokens[i], tokens[i + 1] = '#B', 'A#'
                changed = True
                i += 1
            elif tokens[i] == 'B#' and tokens[i + 1] == '#A':
                tokens[i], tokens[i + 1] = '#A', 'B#'
                changed = True
                i += 1
            elif tokens[i] == 'B#' and tokens[i + 1] == '#B':
                tokens.pop(i)
                tokens.pop(i)
                changed = True
            else:
                i += 1
            
            if changed:
                break  # Restart from the beginning if a change was made
        
        if not changed:
            break  # Exit the loop if no changes were made in the last pass
    
    return tokens

# Initial program
program = ['#A', '#B', '#B', 'B#', '#A', 'A#', '#B', 'B#', '#A', '#A']

# Compute the final state
final_state = compute_program(program)

# Output the final state
print(f"<<<{' '.join(final_state)}>>>")