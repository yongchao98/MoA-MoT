def apply_rules(program):
    i = 0
    while i < len(program) - 1:
        if program[i] == 'A#' and program[i+1] == '#A':
            program.pop(i)
            program.pop(i)
            i = max(i - 1, 0)
        elif program[i] == 'A#' and program[i+1] == '#B':
            program[i], program[i+1] = '#B', 'A#'
            i += 1
        elif program[i] == 'B#' and program[i+1] == '#A':
            program[i], program[i+1] = '#A', 'B#'
            i += 1
        elif program[i] == 'B#' and program[i+1] == '#B':
            program.pop(i)
            program.pop(i)
            i = max(i - 1, 0)
        else:
            i += 1
    return program

initial_program = ['#A', 'A#', '#A', '#B', 'A#', '#B', '#A', 'A#', '#B', 'B#']
final_program = apply_rules(initial_program)
print(' '.join(final_program))