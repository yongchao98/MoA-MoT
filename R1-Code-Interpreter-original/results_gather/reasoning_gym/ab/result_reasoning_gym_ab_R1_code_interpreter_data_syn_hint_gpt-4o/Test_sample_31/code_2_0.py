def apply_rules(program):
    i = 0
    while i < len(program) - 1:
        if program[i] == 'A#' and program[i + 1] == '#A':
            program = program[:i] + program[i + 2:]
            i = max(i - 1, 0)
        elif program[i] == 'A#' and program[i + 1] == '#B':
            program[i], program[i + 1] = '#B', 'A#'
            i += 1
        elif program[i] == 'B#' and program[i + 1] == '#A':
            program[i], program[i + 1] = '#A', 'B#'
            i += 1
        elif program[i] == 'B#' and program[i + 1] == '#B':
            program = program[:i] + program[i + 2:]
            i = max(i - 1, 0)
        else:
            i += 1
    return program

# Initial program
program = ['#B', '#B', '#A', '#A', '#B', 'A#', '#B', 'A#', 'A#', 'B#']

# Apply rules until no more changes can be made
previous_program = None
while program != previous_program:
    previous_program = program[:]
    program = apply_rules(program)

# Print the final state of the program
print(' '.join(program))