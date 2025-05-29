word = "cryptanalytically"
positions = [i + 1 for i, char in enumerate(word) if char == 'y']
count = len(positions)
print(f"Count: {count}, Positions: {positions}")