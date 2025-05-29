def simulate():
    blocks = {
        '[A]': 1, '[B]': 6, '[C]': 6,
        '{A}': 0, '{B}': 0, '{C}': 0,
        '(A)': 0, '(B)': 0, '(C)': 0
    }
    
    def print_state(step, rule):
        state = ", ".join(f"{v} {k}" for k, v in sorted(blocks.items()) if v > 0)
        print(f"After step {step} (Rule {rule}): {state}")
    
    step = 0
    print_state(step, "start")
    
    # First try Rule 3 as much as possible
    while blocks['[B]'] >= 1 and blocks['[C]'] >= 1:
        blocks['[B]'] -= 1
        blocks['[C]'] -= 1
        blocks['{B}'] += 1
        step += 1
        print_state(step, "3")
    
    # Then Rule 4 for remaining [C]s
    while blocks['[C]'] >= 2:
        blocks['[C]'] -= 2
        blocks['{C}'] += 1
        step += 1
        print_state(step, "4")
    
    # Then Rule 2 with remaining [A]
    while blocks['[A]'] >= 1 and blocks['[B]'] >= 1:
        blocks['[A]'] -= 1
        blocks['[B]'] -= 1
        blocks['{C}'] += 1
        step += 1
        print_state(step, "2")
    
    # Then Rule 6 with {B}s
    while blocks['{B}'] >= 2:
        blocks['{B}'] -= 2
        blocks['(C)'] += 1
        step += 1
        print_state(step, "6")
    
    # Finally Rule 5
    while blocks['{A}'] >= 1 and blocks['{C}'] >= 1:
        blocks['{A}'] -= 1
        blocks['{C}'] -= 1
        blocks['(A)'] += 1
        blocks['(B)'] += 1
        step += 1
        print_state(step, "5")
    
    result = ""
    for block, count in sorted(blocks.items()):
        if count > 0:
            result += f"{count} {block}, "
    print("\nFinal state:", result[:-2])

simulate()