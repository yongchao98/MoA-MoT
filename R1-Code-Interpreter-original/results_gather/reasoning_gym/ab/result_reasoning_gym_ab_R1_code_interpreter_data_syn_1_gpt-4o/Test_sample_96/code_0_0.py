def compute_program(tokens):
    i = 0
    while i < len(tokens) - 1:
        pair = tokens[i] + " " + tokens[i + 1]
        if pair == "A# #A":
            tokens.pop(i)
            tokens.pop(i)
            i = max(i - 1, 0)
        elif pair == "A# #B":
            tokens[i], tokens[i + 1] = "#B", "A#"
            i += 1
        elif pair == "B# #A":
            tokens[i], tokens[i + 1] = "#A", "B#"
            i += 1
        elif pair == "B# #B":
            tokens.pop(i)
            tokens.pop(i)
            i = max(i - 1, 0)
        else:
            i += 1
    return tokens

# Initial program
program = ["B#", "#B", "#A", "B#", "#B", "A#", "B#", "#A", "#A", "#A"]

# Compute the final state of the program
final_state = compute_program(program)

# Print the final state
print(final_state)