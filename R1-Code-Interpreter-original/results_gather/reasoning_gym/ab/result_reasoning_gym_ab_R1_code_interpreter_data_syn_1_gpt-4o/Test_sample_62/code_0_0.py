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

def compute_program(initial_tokens):
    tokens = initial_tokens[:]
    while True:
        new_tokens = apply_rules(tokens)
        if new_tokens == tokens:
            break
        tokens = new_tokens
    return tokens

initial_program = ['B#', 'A#', '#A', '#A', '#B', 'A#', '#B', '#A', '#A', '#B']
final_state = compute_program(initial_program)
print(final_state)