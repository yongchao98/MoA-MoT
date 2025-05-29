def apply_rules(program):
    i = 0
    while i < len(program) - 1:
        if program[i] == 'A#' and program[i + 1] == '#A':
            program.pop(i)
            program.pop(i)
        elif program[i] == 'A#' and program[i + 1] == '#B':
            program[i], program[i + 1] = '#B', 'A#'
            i += 1
        elif program[i] == 'B#' and program[i + 1] == '#A':
            program[i], program[i + 1] = '#A', 'B#'
            i += 1
        elif program[i] == 'B#' and program[i + 1] == '#B':
            program.pop(i)
            program.pop(i)
        else:
            i += 1
    return program

# Initial program
program = ['B#', 'B#', 'B#', 'B#', '#A', '#B', '#B', '#A', '#A', '#B']

# Apply rules until no more changes
previous_length = -1
while previous_length != len(program):
    previous_length = len(program)
    program = apply_rules(program)

print(program)