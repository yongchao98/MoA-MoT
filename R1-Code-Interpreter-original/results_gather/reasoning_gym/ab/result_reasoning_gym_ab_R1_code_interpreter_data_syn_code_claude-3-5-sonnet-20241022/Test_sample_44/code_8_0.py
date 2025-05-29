def apply_rules(program):
    tokens = program.split()
    print("Initial:", ' '.join(tokens))
    
    step = 1
    while True:
        changed = False
        
        # First, identify all B# tokens and #A tokens
        b_hash_positions = [i for i, t in enumerate(tokens) if t == "B#"]
        hash_a_positions = [i for i, t in enumerate(tokens) if t == "#A"]
        
        # If we have both B# and #A, they should eventually transform
        if b_hash_positions and hash_a_positions:
            print(f"Step {step}: Found B# at {b_hash_positions} and #A at {hash_a_positions}")
            
        # Process adjacent pairs
        i = 0
        while i < len(tokens) - 1:
            # Check for cancellations first
            if (tokens[i] == "A#" and tokens[i+1] == "#A") or \
               (tokens[i] == "B#" and tokens[i+1] == "#B"):
                print(f"Step {step}: Cancel {tokens[i]} {tokens[i+1]} at {i},{i+1}")
                tokens.pop(i)
                tokens.pop(i)
                print("After:", ' '.join(tokens))
                changed = True
                break
            
            # Check for transformations
            if tokens[i] == "B#" and tokens[i+1] == "#A":
                print(f"Step {step}: Transform B# #A to #A B# at {i},{i+1}")
                tokens[i] = "#A"
                tokens[i+1] = "B#"
                print("After:", ' '.join(tokens))
                changed = True
                break
            
            if tokens[i] == "A#" and tokens[i+1] == "#B":
                print(f"Step {step}: Transform A# #B to #B A# at {i},{i+1}")
                tokens[i] = "#B"
                tokens[i+1] = "A#"
                print("After:", ' '.join(tokens))
                changed = True
                break
            
            i += 1
        
        if not changed:
            # Double check for any B# #A patterns that might be separated
            for i in range(len(tokens)):
                if tokens[i] == "B#":
                    for j in range(i+1, len(tokens)):
                        if tokens[j] == "#A":
                            # We should be able to transform these
                            changed = True
                            break
                    if changed:
                        break
        
        if not changed:
            break
        step += 1
    
    return ' '.join(tokens)

program = "#A #B #B B# #A A# #B B# #A #A"
final_state = apply_rules(program)
print("\nFinal state:", final_state)