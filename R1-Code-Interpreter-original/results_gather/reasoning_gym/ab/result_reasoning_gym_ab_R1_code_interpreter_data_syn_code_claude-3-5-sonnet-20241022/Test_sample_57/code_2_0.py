def apply_rules(program):
    tokens = program.split()
    steps = 1
    
    while True:
        print(f"Step {steps}: {' '.join(tokens)}")
        steps += 1
        
        # Try to find any applicable rule
        found = False
        for i in range(len(tokens) - 1):
            pair = f"{tokens[i]} {tokens[i+1]}"
            
            # Apply rules in order
            if pair == "A# #A":
                tokens[i:i+2] = []  # Remove both tokens
                found = True
                break
            elif pair == "B# #B":
                tokens[i:i+2] = []  # Remove both tokens
                found = True
                break
            elif pair == "A# #B":
                tokens[i:i+2] = ["#B", "A#"]
                found = True
                break
            elif pair == "B# #A":
                tokens[i:i+2] = ["#A", "B#"]
                found = True
                break
        
        if not found:
            break
    
    return " ".join(tokens)

program = "A# #A B# #A A# B# #B #A B# #A"
print(f"Initial: {program}")
print(f"Final: {apply_rules(program)}")