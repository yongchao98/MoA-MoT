def apply_rules(tokens):
    def tokens_to_str(t):
        return ' '.join(t)
    
    print("Initial:", tokens_to_str(tokens))
    
    while True:
        changed = False
        new_tokens = []
        i = 0
        
        while i < len(tokens):
            if i == len(tokens) - 1:
                new_tokens.append(tokens[i])
                i += 1
                continue
                
            pair = tokens[i] + " " + tokens[i+1]
            
            if pair == "A# #A" or pair == "B# #B":
                i += 2
                changed = True
                print("After removal:", tokens_to_str(new_tokens + tokens[i:]))
            elif pair == "A# #B":
                new_tokens.extend(["#B", "A#"])
                i += 2
                changed = True
                print("After swap:", tokens_to_str(new_tokens + tokens[i:]))
            elif pair == "B# #A":
                new_tokens.extend(["#A", "B#"])
                i += 2
                changed = True
                print("After swap:", tokens_to_str(new_tokens + tokens[i:]))
            else:
                new_tokens.append(tokens[i])
                i += 1
        
        if not changed:
            break
        tokens = new_tokens + (tokens[i:] if i < len(tokens) else [])
    
    return tokens_to_str(tokens)

program = "#A #A #B B# A# #A A# A# B# #B".split()
final = apply_rules(program)
print("\nFinal state:", final)