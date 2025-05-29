def apply_rules(program):
    def show_state(tokens, msg=""):
        print(f"{msg}{'Current: ' if msg else 'Initial: '}{' '.join(tokens)}")
    
    tokens = program.split()
    show_state(tokens)
    
    while True:
        i = 0
        changed = False
        while i < len(tokens) - 1:
            # Check if we have a facing ## pattern
            if tokens[i].endswith('#') and tokens[i+1].startswith('#'):
                pair = f"{tokens[i]} {tokens[i+1]}"
                
                if pair in ["B# #B", "A# #A"]:
                    # Delete both tokens
                    tokens = tokens[:i] + tokens[i+2:]
                    changed = True
                    show_state(tokens, f"After {pair} deletion: ")
                    break
                    
                elif pair == "A# #B":
                    tokens = tokens[:i] + ["#B", "A#"] + tokens[i+2:]
                    changed = True
                    show_state(tokens, "After A# #B swap: ")
                    break
                    
                elif pair == "B# #A":
                    tokens = tokens[:i] + ["#A", "B#"] + tokens[i+2:]
                    changed = True
                    show_state(tokens, "After B# #A swap: ")
                    break
            i += 1
            
        if not changed:
            break
    
    return tokens

program = "#B #B #A #A #B A# #B A# A# B#"
final_tokens = apply_rules(program)
print("\nFinal:", " ".join(final_tokens))