def can_apply_any_rule(blocks):
    # Rule 5: {A} + {C} -> (A) + (B)
    if blocks['{A}'] >= 1 and blocks['{C}'] >= 1:
        return True
    # Rule 6: {B} + {B} -> (C)
    if blocks['{B}'] >= 2:
        return True
    return False

blocks = {'{A}': 4, '{B}': 1, '{C}': 0, '(A)': 0, '(B)': 0, '(C)': 0}
if can_apply_any_rule(blocks):
    print("More rules can be applied")
else:
    print("Final state reached")
    print(", ".join(f"{k}: {v}" for k, v in sorted(blocks.items()) if v > 0))