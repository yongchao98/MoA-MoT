word = "counterclassifications"
positions = [i + 1 for i, char in enumerate(word) if char == 's']
count = len(positions)
print(f"Count: {count}, Positions: {positions}")