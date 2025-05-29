def apply_rules(program):
    i = 0
    while i < len(program) - 1:
        if program[i] == 'A#' and program[i+1] == '#A':
            program = program[:i] + program[i+2:]
            i = max(i - 1, 0)
        elif program[i] == 'A#' and program[i+1] == '#B':
            program[i], program[i+1] = program[i+1], program[i]
            i += 1
        elif program[i] == 'B#' and program[i+1] == '#A':
            program[i], program[i+1] = program[i+1], program[i]
            i += 1
        elif program[i] == 'B#' and program[i+1] == '#B':
            program = program[:i] + program[i+2:]
            i = max(i - 1, 0)
        else:
            i += 1
    return program

program = ['A#', '#A', '#B', 'B#', '#B', '#A', '#A', 'A#', '#B', 'B#']
final_state = apply_rules(program)
print(final_state)