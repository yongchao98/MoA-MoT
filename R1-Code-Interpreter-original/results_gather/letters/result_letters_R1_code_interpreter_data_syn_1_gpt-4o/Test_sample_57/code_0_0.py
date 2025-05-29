word = 'nmfghmzgvjinzgmijkdxs'
positions = [i + 1 for i, char in enumerate(word) if char == 'g']
count = len(positions)
print(f"Count: {count}, Positions: {positions}")