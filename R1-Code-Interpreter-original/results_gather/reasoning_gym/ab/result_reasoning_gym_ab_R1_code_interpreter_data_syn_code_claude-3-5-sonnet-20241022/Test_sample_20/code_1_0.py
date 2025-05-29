def apply_rules(program):
    tokens = program.split()
    
    while True:
        changed = False
        i = 0
        while i < len(tokens) - 1:
            pair = tokens[i] + " " + tokens[i+1]
            
            # Apply rules
            if pair in ["A# #A", "B# #B"]:
                tokens.pop(i)
                tokens.pop(i)
                changed = True
                i -= 1
            elif pair == "A# #B":
                tokens[i:i+2] = ["#B", "A#"]
                changed = True
            elif pair == "B# #A":
                tokens[i:i+2] = ["#A", "B#"]
                changed = True
            # Add reverse checks
            elif pair == "#B A#":
                tokens[i:i+2] = ["B#", "#A"]
                changed = True
            elif pair == "#A B#":
                tokens[i:i+2] = ["A#", "#B"]
                changed = True
            i += 1
            
        if not changed:
            break
    
    return " ".join(tokens)

program = "#A A# #B #A #B A# B# B# #B B#"
print(f"Initial: {program}")
print(f"Final: {apply_rules(program)}")