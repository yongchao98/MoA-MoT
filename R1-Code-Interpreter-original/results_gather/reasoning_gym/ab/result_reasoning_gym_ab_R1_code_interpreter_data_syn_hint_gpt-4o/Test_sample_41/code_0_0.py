def apply_rules(program):
    i = 0
    while i < len(program) - 1:
        pair = program[i] + " " + program[i + 1]
        if pair == "A# #A":
            program.pop(i)
            program.pop(i)
        elif pair == "A# #B":
            program[i], program[i + 1] = "#B", "A#"
            i += 1
        elif pair == "B# #A":
            program[i], program[i + 1] = "#A", "B#"
            i += 1
        elif pair == "B# #B":
            program.pop(i)
            program.pop(i)
        else:
            i += 1
    return program

# Initial program
program = ["#A", "A#", "#B", "B#", "#B", "#A", "A#", "#B", "#B", "B#"]

# Apply rules until no more changes can be made
previous_length = -1
while previous_length != len(program):
    previous_length = len(program)
    program = apply_rules(program)

# Print the final state of the program
print(program)