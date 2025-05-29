def apply_method_1(blocks):
    # [A] + [B] + [C] -> {A}
    if blocks['[A]'] >= 1 and blocks['[B]'] >= 1 and blocks['[C]'] >= 1:
        blocks['[A]'] -= 1
        blocks['[B]'] -= 1
        blocks['[C]'] -= 1
        blocks['{A}'] += 1
        return "Method 1: [A]+[B]+[C] -> {A}"
    return None

def apply_method_2(blocks):
    # [A] + [B] -> {C}
    if blocks['[A]'] >= 1 and blocks['[B]'] >= 1:
        blocks['[A]'] -= 1
        blocks['[B]'] -= 1
        blocks['{C}'] += 1
        return "Method 2: [A]+[B] -> {C}"
    return None

def apply_method_3(blocks):
    # [B] + [C] -> {B}
    if blocks['[B]'] >= 1 and blocks['[C]'] >= 1:
        blocks['[B]'] -= 1
        blocks['[C]'] -= 1
        blocks['{B}'] += 1
        return "Method 3: [B]+[C] -> {B}"
    return None

def apply_method_4(blocks):
    # [C] + [C] -> {C}
    if blocks['[C]'] >= 2:
        blocks['[C]'] -= 2
        blocks['{C}'] += 1
        return "Method 4: [C]+[C] -> {C}"
    return None

def apply_method_5(blocks):
    # {A} + {C} -> (A) + (B)
    if blocks['{A}'] >= 1 and blocks['{C}'] >= 1:
        blocks['{A}'] -= 1
        blocks['{C}'] -= 1
        blocks['(A)'] += 1
        blocks['(B)'] += 1
        return "Method 5: {A}+{C} -> (A)+(B)"
    return None

def apply_method_6(blocks):
    # {B} + {B} -> (C)
    if blocks['{B}'] >= 2:
        blocks['{B}'] -= 2
        blocks['(C)'] += 1
        return "Method 6: {B}+{B} -> (C)"
    return None

def simulate():
    blocks = {
        '[A]': 5, '[B]': 4, '[C]': 7,
        '{A}': 0, '{B}': 0, '{C}': 0,
        '(A)': 0, '(B)': 0, '(C)': 0
    }
    
    print("Initial:", blocks)
    step = 1
    
    while True:
        any_change = False
        for method in [apply_method_1, apply_method_2, apply_method_3, 
                      apply_method_4, apply_method_5, apply_method_6]:
            result = method(blocks)
            if result:
                print(f"Step {step}: {result}")
                print("Current state:", {k:v for k,v in blocks.items() if v > 0})
                any_change = True
                step += 1
        
        if not any_change:
            break
    
    print("\nFinal state:")
    result = ""
    for block, count in sorted(blocks.items()):
        if count > 0:
            result += f"{block}:{count} "
    print(result.strip())

simulate()