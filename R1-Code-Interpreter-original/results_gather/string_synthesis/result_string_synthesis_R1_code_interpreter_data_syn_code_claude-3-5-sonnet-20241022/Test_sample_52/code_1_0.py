def synthesize():
    blocks = {
        '[A]': 3, '[B]': 3, '[C]': 5,
        '{A}': 0, '{B}': 0, '{C}': 0,
        '(A)': 0, '(B)': 0, '(C)': 0
    }
    
    while True:
        changes_made = False
        # Try each rule multiple times before moving to next rule
        for rule in [apply_rule1, apply_rule2, apply_rule3, apply_rule4, apply_rule5, apply_rule6]:
            while rule(blocks):
                changes_made = True
        if not changes_made:
            break
    
    result = ""
    for block, count in sorted(blocks.items()):
        if count > 0:
            result += f"{count}{block} "
    print(result.strip())

synthesize()