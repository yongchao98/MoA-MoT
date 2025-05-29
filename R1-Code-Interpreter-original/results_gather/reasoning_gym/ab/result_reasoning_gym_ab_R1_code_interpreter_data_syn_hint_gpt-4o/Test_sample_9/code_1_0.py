def compute_program(tokens):
    changed = True
    while changed:
        changed = False
        i = 0
        while i < len(tokens) - 1:
            pair = tokens[i] + " " + tokens[i + 1]
            if pair == "A# #A":
                tokens = tokens[:i] + tokens[i+2:]
                changed = True
            elif pair == "A# #B":
                tokens[i], tokens[i + 1] = "#B", "A#"
                changed = True
                i += 1
            elif pair == "B# #A":
                tokens[i], tokens[i + 1] = "#A", "B#"
                changed = True
                i += 1
            elif pair == "B# #B":
                tokens = tokens[:i] + tokens[i+2:]
                changed = True
            else:
                i += 1
            if changed:
                break  # Restart from the beginning if any change was made
    return tokens

# Initial program
program = "#B #A B# B# #B #A #A #B A# B#".split()

# Compute the final state
final_state = compute_program(program)

# Print the final state
print(" ".join(final_state))