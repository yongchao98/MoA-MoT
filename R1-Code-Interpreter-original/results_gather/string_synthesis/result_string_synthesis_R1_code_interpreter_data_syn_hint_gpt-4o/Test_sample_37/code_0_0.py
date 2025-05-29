# Final block counts
block_counts = {
    '{A}': 1,
    '(A)': 1,
    '(B)': 1
}

# Format the answer
answer = ''.join(f"{count} {block}" for block, count in block_counts.items())
print(answer)