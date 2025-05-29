def apply_rules(tokens):
    i = 0
    while i < len(tokens) - 1:
        if tokens[i] == 'A#' and tokens[i+1] == '#A':
            tokens.pop(i)
            tokens.pop(i)
        elif tokens[i] == 'A#' and tokens[i+1] == '#B':
            tokens[i], tokens[i+1] = '#B', 'A#'
            i += 1
        elif tokens[i] == 'B#' and tokens[i+1] == '#A':
            tokens[i], tokens[i+1] = '#A', 'B#'
            i += 1
        elif tokens[i] == 'B#' and tokens[i+1] == '#B':
            tokens.pop(i)
            tokens.pop(i)
        else:
            i += 1
    return tokens

# Initial program
program = ['#A', '#A', '#B', 'B#', 'A#', '#A', 'A#', 'A#', 'B#', '#B']

# Apply rules until no more changes
previous_length = -1
while previous_length != len(program):
    previous_length = len(program)
    program = apply_rules(program)

# Output the final state of the program
print(' '.join(program))