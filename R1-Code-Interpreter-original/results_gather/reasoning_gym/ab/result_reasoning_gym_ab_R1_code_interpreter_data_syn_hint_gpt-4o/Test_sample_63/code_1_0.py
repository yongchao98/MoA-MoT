def apply_rules(program):
    changed = True
    while changed:
        changed = False
        i = 0
        while i < len(program) - 1:
            if program[i] == 'A#' and program[i+1] == '#A':
                program = program[:i] + program[i+2:]
                changed = True
            elif program[i] == 'A#' and program[i+1] == '#B':
                program = program[:i] + ['#B', 'A#'] + program[i+2:]
                changed = True
                i += 1
            elif program[i] == 'B#' and program[i+1] == '#A':
                program = program[:i] + ['#A', 'B#'] + program[i+2:]
                changed = True
                i += 1
            elif program[i] == 'B#' and program[i+1] == '#B':
                program = program[:i] + program[i+2:]
                changed = True
            else:
                i += 1
    return program

program = ['B#', '#B', 'A#', 'B#', 'A#', 'A#', '#A', 'B#', 'A#', 'B#']
final_state = apply_rules(program)
print(final_state)