def apply_rules(program):
    i = 0
    while i < len(program) - 1:
        if program[i] == 'A#' and program[i+1] == '#A':
            # Apply A# #A -> nothing
            program.pop(i)
            program.pop(i)
            i = max(i - 1, 0)  # Step back to check for new possible rewrites
        elif program[i] == 'A#' and program[i+1] == '#B':
            # Apply A# #B -> #B A#
            program[i], program[i+1] = '#B', 'A#'
            i += 1
        elif program[i] == 'B#' and program[i+1] == '#A':
            # Apply B# #A -> #A B#
            program[i], program[i+1] = '#A', 'B#'
            i += 1
        elif program[i] == 'B#' and program[i+1] == '#B':
            # Apply B# #B -> nothing
            program.pop(i)
            program.pop(i)
            i = max(i - 1, 0)  # Step back to check for new possible rewrites
        else:
            i += 1
    return program

# Initial program
program = ['#A', 'A#', 'A#', '#A', 'A#', 'A#', 'B#', 'B#', '#A', '#B']

# Apply the rules until no more changes can be made
final_program = apply_rules(program)

# Print the final state of the program
print(final_program)