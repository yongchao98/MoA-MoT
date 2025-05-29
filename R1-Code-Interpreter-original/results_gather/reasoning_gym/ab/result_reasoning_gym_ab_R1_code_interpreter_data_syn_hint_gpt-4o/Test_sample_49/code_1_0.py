def apply_rules(program):
    i = 0
    while i < len(program) - 1:
        pair = program[i] + " " + program[i + 1]
        if pair == "A# #A" or pair == "B# #B":
            # Remove the pair
            program = program[:i] + program[i+2:]
            i = max(i - 1, 0)  # Step back to recheck the previous position
        elif pair == "A# #B":
            # Swap to #B A#
            program = program[:i] + ["#B", "A#"] + program[i+2:]
            i += 1
        elif pair == "B# #A":
            # Swap to #A B#
            program = program[:i] + ["#A", "B#"] + program[i+2:]
            i += 1
        else:
            i += 1
    return program

initial_program = ["#B", "#A", "#B", "#B", "A#", "B#", "B#", "#B", "A#", "#A"]
final_program = apply_rules(initial_program)
print(" ".join(final_program))