def apply_rules(tokens):
    # Return True if any rule was applied, False otherwise
    for i in range(len(tokens) - 1):
        # Check if we have two adjacent tokens
        if i >= len(tokens) - 1:
            break
            
        # Get current pair
        pair = tokens[i:i+2]
        
        # Rule 1: A# #A -> nothing
        if pair == ['A#', '#A']:
            tokens[i:i+2] = []
            return True
            
        # Rule 2: A# #B -> #B A#
        elif pair == ['A#', '#B']:
            tokens[i:i+2] = ['#B', 'A#']
            return True
            
        # Rule 3: B# #A -> #A B#
        elif pair == ['B#', '#A']:
            tokens[i:i+2] = ['#A', 'B#']
            return True
            
        # Rule 4: B# #B -> nothing
        elif pair == ['B#', '#B']:
            tokens[i:i+2] = []
            return True
            
    return False

# Initial program
program = ['B#', '#B', '#A', 'A#', 'A#', 'A#', 'A#', '#A', '#A', 'A#']

# Keep applying rules until no more changes
while True:
    print("Current state:", ' '.join(program))
    if not apply_rules(program):
        break

print("\nFinal state:", ' '.join(program))