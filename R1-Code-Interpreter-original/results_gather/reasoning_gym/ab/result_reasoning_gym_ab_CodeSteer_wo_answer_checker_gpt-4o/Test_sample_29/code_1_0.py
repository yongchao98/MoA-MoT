def compute_program(tokens):
    i = 0
    while i < len(tokens) - 1:
        if tokens[i] == "A#" and tokens[i + 1] == "#A":
            tokens.pop(i)
            tokens.pop(i)
            i = max(i - 1, 0)
        elif tokens[i] == "A#" and tokens[i + 1] == "#B":
            tokens[i], tokens[i + 1] = "#B", "A#"
            i += 1
        elif tokens[i] == "B#" and tokens[i + 1] == "#A":
            tokens[i], tokens[i + 1] = "#A", "B#"
            i += 1
        elif tokens[i] == "B#" and tokens[i + 1] == "#B":
            tokens.pop(i)
            tokens.pop(i)
            i = max(i - 1, 0)
        else:
            i += 1
    return tokens

# Initial program
program = ["#B", "B#", "B#", "#A", "#A", "#A", "A#", "#B", "#A", "#A"]

# Compute the final state
final_state = compute_program(program)

# Output the final state
print(f"<<<{' '.join(final_state)}>>>")