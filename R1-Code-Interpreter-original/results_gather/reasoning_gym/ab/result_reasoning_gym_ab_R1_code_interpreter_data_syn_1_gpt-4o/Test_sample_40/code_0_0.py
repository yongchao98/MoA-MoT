def compute_program(tokens):
    i = 0
    while i < len(tokens) - 1:
        if tokens[i] == 'A#' and tokens[i + 1] == '#A':
            tokens.pop(i)
            tokens.pop(i)
            i = max(i - 1, 0)
        elif tokens[i] == 'A#' and tokens[i + 1] == '#B':
            tokens[i], tokens[i + 1] = '#B', 'A#'
            i += 1
        elif tokens[i] == 'B#' and tokens[i + 1] == '#A':
            tokens[i], tokens[i + 1] = '#A', 'B#'
            i += 1
        elif tokens[i] == 'B#' and tokens[i + 1] == '#B':
            tokens.pop(i)
            tokens.pop(i)
            i = max(i - 1, 0)
        else:
            i += 1
    return tokens

initial_program = ['A#', 'B#', '#A', '#A', '#A', 'B#', '#A', 'B#', '#B', 'B#']
final_state = compute_program(initial_program)
print(final_state)