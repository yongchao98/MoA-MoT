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

# Initial sequence of tokens
tokens = ['B#', '#B', '#A', 'B#', 'B#', 'B#', '#B', 'B#', 'A#', '#A']

# Apply rules until no more changes can be made
previous_length = -1
while previous_length != len(tokens):
    previous_length = len(tokens)
    tokens = apply_rules(tokens)

# Print the final state of the program
print(tokens)