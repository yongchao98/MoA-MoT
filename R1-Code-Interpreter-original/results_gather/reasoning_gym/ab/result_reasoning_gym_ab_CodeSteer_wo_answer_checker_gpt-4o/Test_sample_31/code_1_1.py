def compute_program(tokens):
    i = 0
    while i < len(tokens) - 1:
        if tokens[i] == "A#" and tokens[i + 1] == "#A":
            # Remove both tokens
            tokens.pop(i)
            tokens.pop(i)
            # Move back one step to check the new pair at this position
            i = max(i - 1, 0)
        elif tokens[i] == "A#" and tokens[i + 1] == "#B":
            # Swap them to #B A#
            tokens[i], tokens[i + 1] = "#B", "A#"
            i += 1
        elif tokens[i] == "B#" and tokens[i + 1] == "#A":
            # Swap them to #A B#
            tokens[i], tokens[i + 1] = "#A", "B#"
            i += 1
        elif tokens[i] == "B#" and tokens[i + 1] == "#B":
            # Remove both tokens
            tokens.pop(i)
            tokens.pop(i)
            # Move back one step to check the new pair at this position
            i = max(i - 1, 0)
        else:
            i += 1
    return tokens

# Initial program
program = ["#B", "#B", "#A", "#A", "#B", "A#", "#B", "A#", "A#", "B#"]

# Compute the final state
final_state = compute_program(program)

# Output the final state
print(f"<<<{' '.join(final_state)}>>>")