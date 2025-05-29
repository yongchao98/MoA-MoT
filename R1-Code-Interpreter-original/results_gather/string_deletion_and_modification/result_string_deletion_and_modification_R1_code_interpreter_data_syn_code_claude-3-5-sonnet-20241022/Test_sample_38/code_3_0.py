def find_all_occurrences(s, pattern):
    print(f"Finding all '{pattern}' in '{s}':")
    positions = []
    start = 0
    while True:
        pos = s.find(pattern, start)
        if pos == -1:
            break
        positions.append(pos)
        start = pos + 1
    print(f"Found at positions: {positions}")
    
    # Show context for each occurrence
    for pos in positions:
        before = s[max(0,pos-2):pos]
        after = s[pos+len(pattern):min(len(s),pos+len(pattern)+2)]
        print(f"Position {pos}: ...{before}[{s[pos:pos+len(pattern)]}]{after}...")

print("String analysis for:", "bacbcbacc")
print("Length:", len("bacbcbacc"))
print("\nCharacter by character:")
for i, c in enumerate("bacbcbacc"):
    print(f"Position {i}: {c}")

print("\nSearching for specific patterns:")
patterns = ['ca', 'bca', 'ba', 'cb', 'aa']
for pattern in patterns:
    find_all_occurrences("bacbcbacc", pattern)