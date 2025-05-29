def apply_rules(tokens):
    i = 0
    while i < len(tokens) - 1:
        if tokens[i] == 'A#' and tokens[i + 1] == '#A':
            tokens.pop(i)
            tokens.pop(i)
        elif tokens[i] == 'A#' and tokens[i + 1] == '#B':
            tokens[i], tokens[i + 1] = '#B', 'A#'
            i += 1
        elif tokens[i] == 'B#' and tokens[i + 1] == '#A':
            tokens[i], tokens[i + 1] = '#A', 'B#'
            i += 1
        elif tokens[i] == 'B#' and tokens[i + 1] == '#B':
            tokens.pop(i)
            tokens.pop(i)
        else:
            i += 1
    return tokens

def compute_program(initial_program):
    tokens = initial_program.split()
    previous_length = -1
    while previous_length != len(tokens):
        previous_length = len(tokens)
        tokens = apply_rules(tokens)
    return ' '.join(tokens)

initial_program = "A# A# #A #B A# #B #A #B #B #B"
final_state = compute_program(initial_program)
print(final_state)