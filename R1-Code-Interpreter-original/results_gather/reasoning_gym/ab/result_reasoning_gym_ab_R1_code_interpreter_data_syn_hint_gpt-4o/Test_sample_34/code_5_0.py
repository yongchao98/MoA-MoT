def apply_rules(program):
    changes = True
    while changes:
        changes = False
        i = 0
        while i < len(program) - 1:
            if program[i] == 'A#' and program[i+1] == '#A':
                program = program[:i] + program[i+2:]
                changes = True
            elif program[i] == 'A#' and program[i+1] == '#B':
                program[i], program[i+1] = program[i+1], program[i]
                changes = True
                i += 1
            elif program[i] == 'B#' and program[i+1] == '#A':
                program[i], program[i+1] = program[i+1], program[i]
                changes = True
                i += 1
            elif program[i] == 'B#' and program[i+1] == '#B':
                program = program[:i] + program[i+2:]
                changes = True
            else:
                i += 1
    return program

program = ['A#', '#A', '#B', 'B#', '#B', '#A', '#A', 'A#', '#B', 'B#']
final_state = apply_rules(program)
print(final_state)