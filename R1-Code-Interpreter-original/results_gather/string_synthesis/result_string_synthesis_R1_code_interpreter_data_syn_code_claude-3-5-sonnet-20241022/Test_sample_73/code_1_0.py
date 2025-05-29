def can_apply_rule(blocks, rule_num):
    if rule_num == 1:
        return blocks['[A]'] >= 1 and blocks['[B]'] >= 1 and blocks['[C]'] >= 1
    elif rule_num == 2:
        return blocks['[A]'] >= 1 and blocks['[B]'] >= 1
    elif rule_num == 3:
        return blocks['[B]'] >= 1 and blocks['[C]'] >= 1
    elif rule_num == 4:
        return blocks['[C]'] >= 2
    elif rule_num == 5:
        return blocks['{A}'] >= 1 and blocks['{C}'] >= 1
    elif rule_num == 6:
        return blocks['{B}'] >= 2
    return False

def apply_rule(blocks, rule_num):
    if rule_num == 1:
        blocks['[A]'] -= 1
        blocks['[B]'] -= 1
        blocks['[C]'] -= 1
        blocks['{A}'] += 1
        return "Rule 1: [A] + [B] + [C] → {A}"
    elif rule_num == 2:
        blocks['[A]'] -= 1
        blocks['[B]'] -= 1
        blocks['{C}'] += 1
        return "Rule 2: [A] + [B] → {C}"
    elif rule_num == 3:
        blocks['[B]'] -= 1
        blocks['[C]'] -= 1
        blocks['{B}'] += 1
        return "Rule 3: [B] + [C] → {B}"
    elif rule_num == 4:
        blocks['[C]'] -= 2
        blocks['{C}'] += 1
        return "Rule 4: [C] + [C] → {C}"
    elif rule_num == 5:
        blocks['{A}'] -= 1
        blocks['{C}'] -= 1
        blocks['(A)'] += 1
        blocks['(B)'] += 1
        return "Rule 5: {A} + {C} → (A) + (B)"
    elif rule_num == 6:
        blocks['{B}'] -= 2
        blocks['(C)'] += 1
        return "Rule 6: {B} + {B} → (C)"

def simulate_synthesis():
    blocks = {
        '[A]': 3, '[B]': 2, '[C]': 2,
        '{A}': 0, '{B}': 0, '{C}': 0,
        '(A)': 0, '(B)': 0, '(C)': 0
    }
    
    step = 1
    print(f"Initial: ", end="")
    print(", ".join([f"{v} {k}" for k, v in blocks.items() if v > 0]))
    
    while True:
        changes = False
        for rule in range(1, 7):
            while can_apply_rule(blocks, rule):
                rule_applied = apply_rule(blocks, rule)
                changes = True
                print(f"\nStep {step}: {rule_applied}")
                print("Result: ", end="")
                print(", ".join([f"{v} {k}" for k, v in blocks.items() if v > 0]))
                step += 1
        if not changes:
            break

simulate_synthesis()