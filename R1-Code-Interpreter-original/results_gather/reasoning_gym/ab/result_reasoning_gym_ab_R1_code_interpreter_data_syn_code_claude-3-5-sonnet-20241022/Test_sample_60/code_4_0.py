def apply_rules(program):
    tokens = program.split()
    step = 0
    
    while True:
        print(f"Step {step}: {' '.join(tokens)}")
        changed = False
        i = 0
        while i < len(tokens) - 1:
            current = tokens[i:i+2]
            
            # Rule 1: A# #A -> nothing
            if current == ["A#", "#A"]:
                print(f"Applying rule: A# #A -> nothing at position {i}")
                tokens.pop(i)
                tokens.pop(i)
                changed = True
                break
                
            # Rule 2: B# #B -> nothing
            elif current == ["B#", "#B"]:
                print(f"Applying rule: B# #B -> nothing at position {i}")
                tokens.pop(i)
                tokens.pop(i)
                changed = True
                break
                
            # Rule 3: A# #B -> #B A#
            elif current == ["A#", "#B"]:
                print(f"Applying rule: A# #B -> #B A# at position {i}")
                tokens[i:i+2] = ["#B", "A#"]
                changed = True
                break
                
            # Rule 4: B# #A -> #A B#
            elif current == ["B#", "#A"]:
                print(f"Applying rule: B# #A -> #A B# at position {i}")
                tokens[i:i+2] = ["#A", "B#"]
                changed = True
                break
            
            i += 1
        
        if not changed:
            break
        step += 1
    
    return ' '.join(tokens)

initial = "#B #A #B A# B# A# A# A#"
print("Initial:", initial)
final = apply_rules(initial)
print("\nFinal:", final)