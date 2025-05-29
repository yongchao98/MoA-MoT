def apply_single_rule(tokens):
    for i in range(len(tokens) - 1):
        pair = tokens[i] + " " + tokens[i+1]
        if pair in ["A# #A", "B# #B"]:
            return tokens[:i] + tokens[i+2:], f"Remove pair: {pair}"
        elif pair == "A# #B":
            return tokens[:i] + ["#B", "A#"] + tokens[i+2:], f"Transform: {pair} -> #B A#"
        elif pair == "B# #A":
            return tokens[:i] + ["#A", "B#"] + tokens[i+2:], f"Transform: {pair} -> #A B#"
    return tokens, None

def show_steps(program):
    tokens = program.split()
    step = 1
    print("Initial:", " ".join(tokens))
    
    while True:
        new_tokens, change = apply_single_rule(tokens)
        if change is None:
            break
        tokens = new_tokens
        print(f"Step {step}: {' '.join(tokens)} ({change})")
        step += 1
    
    print("Final:", " ".join(tokens))

program = "#A B# #B A# A# #B B# B# #A #A"
show_steps(program)