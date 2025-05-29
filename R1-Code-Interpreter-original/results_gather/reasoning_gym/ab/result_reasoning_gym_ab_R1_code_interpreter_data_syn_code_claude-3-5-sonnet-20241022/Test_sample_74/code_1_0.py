def apply_rules(tokens):
    tokens = tokens.split()
    step = 1
    
    while True:
        changed = False
        i = 0
        print(f"\nStep {step}: {' '.join(tokens)}")
        while i < len(tokens) - 1:
            pair = tokens[i] + " " + tokens[i+1]
            print(f"  Checking pair: {pair}")
            
            if pair == "A# #A" or pair == "B# #B":
                print(f"    Deleting pair: {pair}")
                tokens.pop(i)
                tokens.pop(i)
                changed = True
                break
            elif pair == "A# #B":
                print(f"    Transforming: {pair} -> #B A#")
                tokens[i:i+2] = ["#B", "A#"]
                changed = True
                break
            elif pair == "B# #A":
                print(f"    Transforming: {pair} -> #A B#")
                tokens[i:i+2] = ["#A", "B#"]
                changed = True
                break
            i += 1
        
        if not changed:
            break
        step += 1
    
    return " ".join(tokens)

program = "#B A# #B A# #A A# B# #B #A #B"
print(f"Initial: {program}")
result = apply_rules(program)
print(f"\nFinal: {result}")