word = "parallelogrammatical"
positions = [i + 1 for i, char in enumerate(word) if char == 'a']
count = len(positions)
print(f"Count: {count}, Positions: {positions}")