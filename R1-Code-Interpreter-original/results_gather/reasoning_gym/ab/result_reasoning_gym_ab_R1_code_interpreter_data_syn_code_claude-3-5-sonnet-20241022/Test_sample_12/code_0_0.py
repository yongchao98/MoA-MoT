def can_rewrite(token1, token2):
    # Returns the result of rewriting two tokens, or None if no rule applies
    if token1 == "A#" and token2 == "#A":
        return []  # disappears
    elif token1 == "A#" and token2 == "#B":
        return ["#B", "A#"]
    elif token1 == "B#" and token2 == "#A":
        return ["#A", "B#"]
    elif token1 == "B#" and token2 == "#B":
        return []  # disappears
    return None

def rewrite_program(program):
    tokens = program.split()
    changed = True
    
    while changed:
        changed = False
        new_tokens = []
        i = 0
        while i < len(tokens):
            if i == len(tokens) - 1:
                new_tokens.append(tokens[i])
                i += 1
            else:
                result = can_rewrite(tokens[i], tokens[i + 1])
                if result is None:
                    new_tokens.append(tokens[i])
                    i += 1
                else:
                    new_tokens.extend(result)
                    i += 2
                    changed = True
        tokens = new_tokens
    
    return " ".join(tokens)

# Initial program
program = "#B A# B# #B B# #A A# B# A# A#"
final_state = rewrite_program(program)
print(final_state)