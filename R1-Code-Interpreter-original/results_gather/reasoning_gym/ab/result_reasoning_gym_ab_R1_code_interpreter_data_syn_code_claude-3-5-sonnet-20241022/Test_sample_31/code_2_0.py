def apply_rules(program):
    def show(tokens, step=""):
        print(f"{step}: {' '.join(tokens)}")
    
    tokens = program.split()
    show(tokens, "Initial")
    
    while True:
        i = 0
        changed = False
        
        # Convert to list for easier manipulation
        while i < len(tokens) - 1:
            current = tokens[i]
            next_token = tokens[i + 1]
            
            # Check if current pair ends with # and starts with #
            if current.endswith('#') and next_token.startswith('#'):
                if (current == "A#" and next_token == "#A") or (current == "B#" and next_token == "#B"):
                    # Remove both tokens
                    tokens = tokens[:i] + tokens[i+2:]
                    changed = True
                    show(tokens, "After deletion")
                    break
                elif current == "A#" and next_token == "#B":
                    # Replace "A# #B" with "#B A#"
                    tokens = tokens[:i] + ["#B", "A#"] + tokens[i+2:]
                    changed = True
                    show(tokens, "After A#-#B swap")
                    break
                elif current == "B#" and next_token == "#A":
                    # Replace "B# #A" with "#A B#"
                    tokens = tokens[:i] + ["#A", "B#"] + tokens[i+2:]
                    changed = True
                    show(tokens, "After B#-#A swap")
                    break
            i += 1
        
        if not changed:
            break
    
    return tokens

program = "#B #B #A #A #B A# #B A# A# B#"
final_tokens = apply_rules(program)
print("\nFinal:", " ".join(final_tokens))