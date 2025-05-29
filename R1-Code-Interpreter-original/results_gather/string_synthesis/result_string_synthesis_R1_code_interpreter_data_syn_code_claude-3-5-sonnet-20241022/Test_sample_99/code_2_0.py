def simulate_final(blocks):
    # Try Method 3: [B] + [C] -> {B}
    if blocks['[B]'] >= 1 and blocks['[C]'] >= 1:
        blocks['[B]'] -= 1
        blocks['[C]'] -= 1
        blocks['{B}'] += 1
        return True
    return False

blocks = {
    '[A]': 0, '[B]': 2, '[C]': 1,
    '{A}': 0, '{B}': 0, '{C}': 5,
    '(A)': 0, '(B)': 0, '(C)': 0
}

changed = True
while changed:
    changed = simulate_final(blocks)

result = ""
for block, count in sorted(blocks.items()):
    if count > 0:
        result += f"{count} {block}, "
print(result.rstrip(", "))