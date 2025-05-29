def can_process(tokens, i):
    if i >= len(tokens) - 1:
        return False
    t1, t2 = tokens[i], tokens[i+1]
    return (t1 == "A#" and t2 == "#A") or \
           (t1 == "A#" and t2 == "#B") or \
           (t1 == "B#" and t2 == "#A") or \
           (t1 == "B#" and t2 == "#B")

def process_pair(tokens, i):
    t1, t2 = tokens[i], tokens[i+1]
    if t1 == "A#" and t2 == "#A":
        return tokens[:i] + tokens[i+2:]
    elif t1 == "A#" and t2 == "#B":
        return tokens[:i] + ["#B", "A#"] + tokens[i+2:]
    elif t1 == "B#" and t2 == "#A":
        return tokens[:i] + ["#A", "B#"] + tokens[i+2:]
    elif t1 == "B#" and t2 == "#B":
        return tokens[:i] + tokens[i+2:]
    return tokens

def process_program(program):
    tokens = program.split()
    while True:
        found_change = False
        for i in range(len(tokens)):
            if can_process(tokens, i):
                tokens = process_pair(tokens, i)
                found_change = True
                break
        if not found_change:
            break
    return " ".join(tokens)

program = "#A B# #A B# #A #B #A A# #A #A"
result = process_program(program)
print(result)