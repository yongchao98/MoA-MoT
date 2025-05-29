def apply_rules(program):
    changed = True
    while changed:
        changed = False
        i = 0
        while i < len(program) - 1:
            if program[i] == 'A#' and program[i+1] == '#A':
                program.pop(i)
                program.pop(i)
                changed = True
            elif program[i] == 'A#' and program[i+1] == '#B':
                program[i], program[i+1] = '#B', 'A#'
                changed = True
                i += 1
            elif program[i] == 'B#' and program[i+1] == '#A':
                program[i], program[i+1] = '#A', 'B#'
                changed = True
                i += 1
            elif program[i] == 'B#' and program[i+1] == '#B':
                program.pop(i)
                program.pop(i)
                changed = True
            else:
                i += 1
    return program

program = ['A#', 'B#', 'B#', '#A', 'A#', 'B#', '#A', 'A#', '#A', '#B']
final_state = apply_rules(program)
print(final_state)