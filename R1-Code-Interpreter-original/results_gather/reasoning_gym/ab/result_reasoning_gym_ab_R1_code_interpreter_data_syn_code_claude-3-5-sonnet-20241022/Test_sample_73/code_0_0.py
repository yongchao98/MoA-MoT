def process_program(program):
    # Split into tokens
    tokens = program.split()
    
    # Keep processing until no changes can be made
    while True:
        changed = False
        new_tokens = []
        i = 0
        
        while i < len(tokens):
            # If this is the last token, just add it
            if i == len(tokens) - 1:
                new_tokens.append(tokens[i])
                i += 1
                continue
                
            # Get current pair of tokens
            current = tokens[i]
            next_token = tokens[i + 1]
            
            # Apply rules
            if current == "A#" and next_token == "#A":
                # Both disappear
                changed = True
                i += 2
            elif current == "A#" and next_token == "#B":
                # Swap positions
                new_tokens.extend(["#B", "A#"])
                changed = True
                i += 2
            elif current == "B#" and next_token == "#A":
                # Swap positions
                new_tokens.extend(["#A", "B#"])
                changed = True
                i += 2
            elif current == "B#" and next_token == "#B":
                # Both disappear
                changed = True
                i += 2
            else:
                # No rule applies, keep first token
                new_tokens.append(current)
                i += 1
        
        # If no changes were made, we're done
        if not changed:
            return " ".join(new_tokens)
        
        tokens = new_tokens

# Initial program
program = "#A A# A# #A #A #A B# #A B# B#"
print(process_program(program))