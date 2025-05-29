def apply_rules_until_stable(program):
    while True:
        new_program = []
        i = 0
        changed = False
        while i < len(program):
            if i < len(program) - 1:
                if program[i] == 'A#' and program[i+1] == '#A':
                    i += 2
                    changed = True
                elif program[i] == 'A#' and program[i+1] == '#B':
                    new_program.extend(['#B', 'A#'])
                    i += 2
                    changed = True
                elif program[i] == 'B#' and program[i+1] == '#A':
                    new_program.extend(['#A', 'B#'])
                    i += 2
                    changed = True
                elif program[i] == 'B#' and program[i+1] == '#B':
                    i += 2
                    changed = True
                else:
                    new_program.append(program[i])
                    i += 1
            else:
                new_program.append(program[i])
                i += 1
        if not changed:
            break
        program = new_program
    return program

# Initial program
program = ['#A', 'A#', '#B', '#A', '#B', 'A#', 'B#', 'B#', '#B', 'B#']

# Apply rules until no more changes can be made
final_program = apply_rules_until_stable(program)

# Print the final state of the program
print(final_program)