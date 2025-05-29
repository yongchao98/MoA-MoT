def apply_rules(program):
    tokens = program.split()
    print("Starting state:", " ".join(tokens))
    
    while True:
        changed = False
        i = 0
        while i < len(tokens) - 1:
            current = tokens[i]
            next_token = tokens[i + 1]
            
            # Print current pair being checked
            print(f"Checking pair: {current} {next_token}")
            
            if current == "A#" and next_token == "#A":
                tokens.pop(i)
                tokens.pop(i)
                changed = True
                print("Applied rule: A# #A -> nothing")
                break
            elif current == "A#" and next_token == "#B":
                tokens[i] = "#B"
                tokens[i + 1] = "A#"
                changed = True
                print("Applied rule: A# #B -> #B A#")
                break
            elif current == "B#" and next_token == "#A":
                tokens[i] = "#A"
                tokens[i + 1] = "B#"
                changed = True
                print("Applied rule: B# #A -> #A B#")
                break
            elif current == "B#" and next_token == "#B":
                tokens.pop(i)
                tokens.pop(i)
                changed = True
                print("Applied rule: B# #B -> nothing")
                break
            
            i += 1
        
        if not changed:
            break
        print("Current state:", " ".join(tokens))
    
    return " ".join(tokens)

program = "#B #A #A B#"
final_state = apply_rules(program)
print("\nFinal state:", final_state)