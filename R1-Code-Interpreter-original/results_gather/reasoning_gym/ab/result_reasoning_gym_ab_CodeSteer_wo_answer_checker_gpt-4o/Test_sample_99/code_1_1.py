def compute_program(tokens):
    result = []
    i = 0
    while i < len(tokens):
        if i < len(tokens) - 1:
            if tokens[i] == 'A#' and tokens[i + 1] == '#A':
                # Apply A# #A -> nothing
                i += 2
            elif tokens[i] == 'A#' and tokens[i + 1] == '#B':
                # Apply A# #B -> #B A#
                result.append('#B')
                result.append('A#')
                i += 2
            elif tokens[i] == 'B#' and tokens[i + 1] == '#A':
                # Apply B# #A -> #A B#
                result.append('#A')
                result.append('B#')
                i += 2
            elif tokens[i] == 'B#' and tokens[i + 1] == '#B':
                # Apply B# #B -> nothing
                i += 2
            else:
                result.append(tokens[i])
                i += 1
        else:
            result.append(tokens[i])
            i += 1

    return result

# Initial program
program = ['#A', '#A', 'A#', '#B', '#A', 'A#', '#B', 'A#', '#A', 'A#']

# Compute the final state
final_state = compute_program(program)

# Output the final state
print(f"<<<{' '.join(final_state)}>>>")