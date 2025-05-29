word = "nondeterministically"
positions = [i + 1 for i, char in enumerate(word) if char == 'm']
count = len(positions)
print(f"<<<Count: {count}, Positions: {positions}>>>")