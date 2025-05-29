# Final count of blocks
blocks = {
    '[A]': 4,
    '[B]': 0,
    '[C]': 0,
    '{A}': 1,
    '{B}': 0,
    '{C}': 0,
    '(A)': 0,
    '(B)': 0,
    '(C)': 0
}

# Format the answer
answer = ''.join(f"{count} {block}" for block, count in blocks.items() if count > 0)
print(answer)