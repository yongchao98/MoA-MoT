def compute_program(tokens):
    i = 0
    while i < len(tokens) - 1:
        if tokens[i] == 'A#' and tokens[i + 1] == '#A':
            tokens = tokens[:i] + tokens[i + 2:]
            i = 0
        elif tokens[i] == 'A#' and tokens[i + 1] == '#B':
            tokens = tokens[:i] + ['#B', 'A#'] + tokens[i + 2:]
            i = 0
        elif tokens[i] == 'B#' and tokens[i + 1] == '#A':
            tokens = tokens[:i] + ['#A', 'B#'] + tokens[i + 2:]
            i = 0
        elif tokens[i] == 'B#' and tokens[i + 1] == '#B':
            tokens = tokens[:i] + tokens[i + 2:]
            i = 0
        else:
            i += 1
    return tokens

initial_program = ['#B', '#B', 'A#', 'A#', '#A', '#B', '#A', 'A#', 'B#', '#A']
final_state = compute_program(initial_program)
print(final_state)