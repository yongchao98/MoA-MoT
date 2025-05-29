def can_reduce(tokens, i):
    if i + 1 >= len(tokens):
        return False
    pair = tokens[i] + " " + tokens[i+1]
    return pair in ["A# #A", "A# #B", "B# #A", "B# #B"]

def apply_rule(tokens, i):
    pair = tokens[i] + " " + tokens[i+1]
    if pair == "A# #A" or pair == "B# #B":
        return tokens[:i] + tokens[i+2:]
    elif pair == "A# #B":
        return tokens[:i] + ["#B", "A#"] + tokens[i+2:]
    elif pair == "B# #A":
        return tokens[:i] + ["#A", "B#"] + tokens[i+2:]
    return tokens

def reduce_program(program):
    tokens = program.split()
    changed = True
    while changed:
        changed = False
        i = 0
        while i < len(tokens):
            if can_reduce(tokens, i):
                tokens = apply_rule(tokens, i)
                changed = True
                # Print intermediate state
                print("->", " ".join(tokens))
            else:
                i += 1
    return " ".join(tokens)

initial = "B# A# B# #B #A #A B# B# B# A#"
print("Initial:", initial)
result = reduce_program(initial)
print("\nFinal:", result)