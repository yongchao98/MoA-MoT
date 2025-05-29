def apply_rules(program):
    i = 0
    while i < len(program) - 1:
        if program[i] == 'A#' and program[i+1] == '#A':
            program = program[:i] + program[i+2:]
            i = 0  # Restart from the beginning
        elif program[i] == 'A#' and program[i+1] == '#B':
            program = program[:i] + ['#B', 'A#'] + program[i+2:]
            i = 0  # Restart from the beginning
        elif program[i] == 'B#' and program[i+1] == '#A':
            program = program[:i] + ['#A', 'B#'] + program[i+2:]
            i = 0  # Restart from the beginning
        elif program[i] == 'B#' and program[i+1] == '#B':
            program = program[:i] + program[i+2:]
            i = 0  # Restart from the beginning
        else:
            i += 1
    return program

# Initial program
program = ['B#', '#B', '#A', 'B#', 'B#', '#A', 'B#', 'A#', '#B', '#B']

# Apply rules until no more changes can be made
final_state = apply_rules(program)

# Print the final state of the program
print(final_state)