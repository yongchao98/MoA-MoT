def simulate():
    blocks = {
        '[A]': 1, '[B]': 6, '[C]': 6,
        '{A}': 0, '{B}': 0, '{C}': 0,
        '(A)': 0, '(B)': 0, '(C)': 0
    }
    
    def print_state(step):
        state = ", ".join(f"{v} {k}" for k, v in blocks.items() if v > 0)
        print(f"After step {step}: {state}")
    
    def try_all_rules():
        # Rule 4: [C] + [C] -> {C}
        while blocks['[C]'] >= 2:
            blocks['[C]'] -= 2
            blocks['{C}'] += 1
            return True
            
        # Rule 2: [A] + [B] -> {C}
        if blocks['[A]'] >= 1 and blocks['[B]'] >= 1:
            blocks['[A]'] -= 1
            blocks['[B]'] -= 1
            blocks['{C}'] += 1
            return True
            
        # Rule 3: [B] + [C] -> {B}
        if blocks['[B]'] >= 1 and blocks['[C]'] >= 1:
            blocks['[B]'] -= 1
            blocks['[C]'] -= 1
            blocks['{B}'] += 1
            return True
            
        # Rule 1: [A] + [B] + [C] -> {A}
        if blocks['[A]'] >= 1 and blocks['[B]'] >= 1 and blocks['[C]'] >= 1:
            blocks['[A]'] -= 1
            blocks['[B]'] -= 1
            blocks['[C]'] -= 1
            blocks['{A}'] += 1
            return True
            
        # Rule 5: {A} + {C} -> (A) + (B)
        if blocks['{A}'] >= 1 and blocks['{C}'] >= 1:
            blocks['{A}'] -= 1
            blocks['{C}'] -= 1
            blocks['(A)'] += 1
            blocks['(B)'] += 1
            return True
            
        # Rule 6: {B} + {B} -> (C)
        if blocks['{B}'] >= 2:
            blocks['{B}'] -= 2
            blocks['(C)'] += 1
            return True
            
        return False

    step = 0
    print_state(step)
    
    while try_all_rules():
        step += 1
        print_state(step)
    
    result = ""
    for block, count in sorted(blocks.items()):
        if count > 0:
            result += f"{count} {block}, "
    print("\nFinal state:", result[:-2])

simulate()