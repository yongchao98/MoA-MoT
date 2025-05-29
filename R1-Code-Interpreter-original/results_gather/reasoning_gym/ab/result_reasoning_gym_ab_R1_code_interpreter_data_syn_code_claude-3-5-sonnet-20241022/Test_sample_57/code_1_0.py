def apply_rules(program):
    tokens = program.split()
    steps = 1
    
    while True:
        changed = False
        i = 0
        while i < len(tokens) - 1:
            pair = f"{tokens[i]} {tokens[i+1]}"
            
            # Print each step for verification
            print(f"Step {steps}: {' '.join(tokens)}")
            steps += 1
            
            if pair == "A# #A" or pair == "B# #B":
                tokens.pop(i)
                tokens.pop(i)
                changed = True
                break
            elif pair == "A# #B":
                tokens[i:i+2] = ["#B", "A#"]
                changed = True
                break
            elif pair == "B# #A":
                tokens[i:i+2] = ["#A", "B#"]
                changed = True
                break
            i += 1
            
        if not changed:
            break
    
    return " ".join(tokens)

program = "A# #A B# #A A# B# #B #A B# #A"
print(f"Initial: {program}")
print(f"Final: {apply_rules(program)}")