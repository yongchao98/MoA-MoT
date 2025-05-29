def apply_rules(program):
    i = 0
    while i < len(program) - 1:
        if program[i] == 'A#' and program[i+1] == '#A':
            program = program[:i] + program[i+2:]
            i = max(i - 1, 0)
        elif program[i] == 'A#' and program[i+1] == '#B':
            program = program[:i] + ['#B', 'A#'] + program[i+2:]
            i += 1
        elif program[i] == 'B#' and program[i+1] == '#A':
            program = program[:i] + ['#A', 'B#'] + program[i+2:]
            i += 1
        elif program[i] == 'B#' and program[i+1] == '#B':
            program = program[:i] + program[i+2:]
            i = max(i - 1, 0)
        else:
            i += 1
    return program

# Initial program
program = ['B#', '#B', '#A', 'B#', '#B', 'A#', 'B#', '#A', '#A', '#A']

# Apply rules until no more changes
previous_program = None
while program != previous_program:
    previous_program = program
    program = apply_rules(program)

# Print the final state of the program
print(program)