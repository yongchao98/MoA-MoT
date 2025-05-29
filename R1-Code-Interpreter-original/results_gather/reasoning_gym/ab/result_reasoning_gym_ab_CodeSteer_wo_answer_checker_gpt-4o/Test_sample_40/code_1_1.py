def compute_program(tokens):
    changed = True
    while changed:
        changed = False
        i = 0
        while i < len(tokens) - 1:
            if tokens[i] == "A#" and tokens[i + 1] == "#A":
                # A# #A becomes nothing
                tokens.pop(i)
                tokens.pop(i)
                changed = True
            elif tokens[i] == "B#" and tokens[i + 1] == "#B":
                # B# #B becomes nothing
                tokens.pop(i)
                tokens.pop(i)
                changed = True
            elif tokens[i] == "A#" and tokens[i + 1] == "#B":
                # A# #B becomes #B A#
                tokens[i], tokens[i + 1] = "#B", "A#"
                changed = True
                i += 1
            elif tokens[i] == "B#" and tokens[i + 1] == "#A":
                # B# #A becomes #A B#
                tokens[i], tokens[i + 1] = "#A", "B#"
                changed = True
                i += 1
            else:
                i += 1
            # If a change was made, reset i to recheck from the start
            if changed:
                i = max(i - 1, 0)
    return tokens

# Initial program
program = ["A#", "B#", "#A", "#A", "#A", "B#", "#A", "B#", "#B", "B#"]

# Compute the final state
final_state = compute_program(program)

# Output the final state
print(f"<<<{' '.join(final_state)}>>>")