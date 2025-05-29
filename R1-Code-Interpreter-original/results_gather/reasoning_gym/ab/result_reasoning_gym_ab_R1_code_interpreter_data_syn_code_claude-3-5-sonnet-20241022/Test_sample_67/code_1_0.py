def apply_rules(tokens):
    # Return True if any rule was applied, False otherwise
    
    # Helper function to convert list to string for easier comparison
    def list_to_str(lst):
        return ' '.join(lst)
    
    # Check each pair of adjacent tokens
    i = 0
    while i < len(tokens) - 1:
        pair = tokens[i:i+2]
        
        # Rule: B# #B -> nothing
        if pair == ['B#', '#B']:
            tokens[i:i+2] = []
            return True
            
        # Rule: A# #A -> nothing
        elif pair == ['A#', '#A']:
            tokens[i:i+2] = []
            return True
            
        # Rule: A# #B -> #B A#
        elif pair == ['A#', '#B']:
            tokens[i:i+2] = ['#B', 'A#']
            return True
            
        # Rule: B# #A -> #A B#
        elif pair == ['B#', '#A']:
            tokens[i:i+2] = ['#A', 'B#']
            return True
        
        i += 1
    
    return False

# Initial program
program = ['#B', '#B', 'B#', 'B#', 'B#', 'B#', 'B#', 'A#', '#A', '#B']

# Keep applying rules until no more changes
steps = 0
while True:
    if not apply_rules(program):
        break
    steps += 1
    print(f"Step {steps}: {' '.join(program)}")

print("Final:", ' '.join(program))