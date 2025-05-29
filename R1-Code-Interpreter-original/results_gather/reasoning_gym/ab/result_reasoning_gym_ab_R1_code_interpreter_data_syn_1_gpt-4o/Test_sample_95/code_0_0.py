def apply_rules(tokens):
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

# Initial sequence
tokens = ["A#", "A#", "A#", "A#", "#A", "B#", "#B", "B#", "#A", "#A"]

# Apply rules until no more changes
previous_length = -1
while previous_length != len(tokens):
    previous_length = len(tokens)
    tokens = apply_rules(tokens)

# Print the final state of the program
print(" ".join(tokens))