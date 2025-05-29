def apply_rules(program):
    i = 0
    while i < len(program) - 1:
        if program[i] == 'A#' and program[i + 1] == '#A':
            program = program[:i] + program[i + 2:]
            i = 0  # Restart from the beginning
        elif program[i] == 'A#' and program[i + 1] == '#B':
            program = program[:i] + ['#B', 'A#'] + program[i + 2:]
            i = 0  # Restart from the beginning
        elif program[i] == 'B#' and program[i + 1] == '#A':
            program = program[:i] + ['#A', 'B#'] + program[i + 2:]
            i = 0  # Restart from the beginning
        elif program[i] == 'B#' and program[i + 1] == '#B':
            program = program[:i] + program[i + 2:]
            i = 0  # Restart from the beginning
        else:
            i += 1
    return program

# Initial program
program = ['#B', 'B#', 'B#', '#A', '#A', '#A', 'A#', '#B', '#A', '#A']

# Apply the rules
final_state = apply_rules(program)

# Print the final state
print(final_state)