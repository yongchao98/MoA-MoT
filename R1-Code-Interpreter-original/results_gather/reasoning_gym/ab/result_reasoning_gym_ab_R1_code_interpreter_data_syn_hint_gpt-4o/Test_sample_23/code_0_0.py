def apply_rules(program):
    i = 0
    while i < len(program) - 1:
        if program[i] == 'A#' and program[i+1] == '#A':
            program.pop(i)
            program.pop(i)
        elif program[i] == 'A#' and program[i+1] == '#B':
            program[i], program[i+1] = '#B', 'A#'
            i += 1
        elif program[i] == 'B#' and program[i+1] == '#A':
            program[i], program[i+1] = '#A', 'B#'
            i += 1
        elif program[i] == 'B#' and program[i+1] == '#B':
            program.pop(i)
            program.pop(i)
        else:
            i += 1
    return program

program = ['A#', 'B#', 'B#', '#A', 'A#', 'B#', '#A', 'A#', '#A', '#B']
final_state = apply_rules(program)
print(final_state)