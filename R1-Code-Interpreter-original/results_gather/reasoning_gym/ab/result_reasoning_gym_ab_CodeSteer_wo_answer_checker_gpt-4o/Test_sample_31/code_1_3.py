def compute_program(tokens):
    stack = []
    for token in tokens:
        if stack:
            top = stack[-1]
            if top == "A#" and token == "#A":
                stack.pop()  # Remove both tokens
            elif top == "B#" and token == "#B":
                stack.pop()  # Remove both tokens
            elif top == "A#" and token == "#B":
                stack[-1] = "#B"  # Swap them to #B A#
                stack.append("A#")
            elif top == "B#" and token == "#A":
                stack[-1] = "#A"  # Swap them to #A B#
                stack.append("B#")
            else:
                stack.append(token)
        else:
            stack.append(token)
    return stack

# Initial program
program = ["#B", "#B", "#A", "#A", "#B", "A#", "#B", "A#", "A#", "B#"]

# Compute the final state
final_state = compute_program(program)

# Output the final state
print(f"<<<{' '.join(final_state)}>>>")