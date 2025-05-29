def find_applicable_rule(tokens, pos):
    if pos >= len(tokens) - 1:
        return None
    
    t1, t2 = tokens[pos], tokens[pos+1]
    
    # Check if the '#' symbols face each other
    if t1.startswith('#') and t2.endswith('#'):
        return None
    if t1.endswith('#') and t2.startswith('#'):
        pair = [t1, t2]
        if pair == ['A#', '#A'] or pair == ['B#', '#B']:
            return 'eliminate'
        elif pair == ['A#', '#B']:
            return 'swap_AB'
        elif pair == ['B#', '#A']:
            return 'swap_BA'
    return None

def apply_rule(tokens, pos, rule):
    if rule == 'eliminate':
        return tokens[:pos] + tokens[pos+2:]
    elif rule == 'swap_AB':
        return tokens[:pos] + ['#B', 'A#'] + tokens[pos+2:]
    elif rule == 'swap_BA':
        return tokens[:pos] + ['#A', 'B#'] + tokens[pos+2:]
    return tokens

def analyze_sequence(tokens):
    print(f"Current sequence: {' '.join(tokens)}")
    for i in range(len(tokens)-1):
        rule = find_applicable_rule(tokens, i)
        if rule:
            print(f"  Position {i}: Can apply rule {rule} to {tokens[i]} {tokens[i+1]}")
            return i, rule
    return -1, None

# Initial sequence
tokens = ['#B', '#A', '#A', 'B#', '#B', '#B', '#B', '#A', '#B', 'B#']

step = 0
print(f"\nStep {step}:")
pos, rule = analyze_sequence(tokens)

while rule:
    step += 1
    tokens = apply_rule(tokens, pos, rule)
    print(f"\nStep {step}:")
    pos, rule = analyze_sequence(tokens)

print("\nFinal state reached!")