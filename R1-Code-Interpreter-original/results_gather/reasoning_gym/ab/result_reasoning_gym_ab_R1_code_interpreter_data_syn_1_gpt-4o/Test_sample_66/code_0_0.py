def apply_rules(program):
    i = 0
    while i < len(program) - 1:
        if program[i] == 'A#' and program[i + 1] == '#A':
            # A# #A becomes nothing
            program = program[:i] + program[i+2:]
            i = max(i - 1, 0)  # Step back to check new neighbors
        elif program[i] == 'A#' and program[i + 1] == '#B':
            # A# #B becomes #B A#
            program[i], program[i + 1] = '#B', 'A#'
            i += 1
        elif program[i] == 'B#' and program[i + 1] == '#A':
            # B# #A becomes #A B#
            program[i], program[i + 1] = '#A', 'B#'
            i += 1
        elif program[i] == 'B#' and program[i + 1] == '#B':
            # B# #B becomes nothing
            program = program[:i] + program[i+2:]
            i = max(i - 1, 0)  # Step back to check new neighbors
        else:
            i += 1
    return program

# Initial program
program = ['#A', 'A#', '#A', '#B', 'A#', '#B', '#A', 'A#', '#B', 'B#']

# Apply rules until no more changes
final_state = apply_rules(program)

# Print the final state
print(final_state)