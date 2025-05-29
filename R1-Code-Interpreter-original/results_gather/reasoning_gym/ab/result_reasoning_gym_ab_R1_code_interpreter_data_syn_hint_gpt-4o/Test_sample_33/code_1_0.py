def apply_rules(program):
    changed = True
    while changed:
        changed = False
        i = 0
        while i < len(program) - 1:
            if program[i] == 'A#' and program[i + 1] == '#A':
                program = program[:i] + program[i + 2:]
                changed = True
            elif program[i] == 'A#' and program[i + 1] == '#B':
                program = program[:i] + ['#B', 'A#'] + program[i + 2:]
                changed = True
            elif program[i] == 'B#' and program[i + 1] == '#A':
                program = program[:i] + ['#A', 'B#'] + program[i + 2:]
                changed = True
            elif program[i] == 'B#' and program[i + 1] == '#B':
                program = program[:i] + program[i + 2:]
                changed = True
            if changed:
                break
            i += 1
    return program

initial_program = ['#B', '#A', '#A', 'B#', '#B', '#B', '#B', '#A', '#B', 'B#']
final_program = apply_rules(initial_program)
print(final_program)