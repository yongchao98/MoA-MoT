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
                i += 1
                changed = True
            elif program[i] == 'B#' and program[i+1] == '#A':
                program[i], program[i+1] = '#A', 'B#'
                i += 1
                changed = True
            elif program[i] == 'B#' and program[i+1] == '#B':
                program.pop(i)
                program.pop(i)
                changed = True
            else:
                i += 1
    return program

program = ['B#', '#A', '#B', 'A#', 'B#', 'B#', '#A', '#B', 'B#', '#A']
final_state = apply_rules(program)
print(final_state)