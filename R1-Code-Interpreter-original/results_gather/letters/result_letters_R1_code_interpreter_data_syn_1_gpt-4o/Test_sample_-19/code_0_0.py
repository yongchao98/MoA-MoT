word = "stereoroentgenography"
positions = [i + 1 for i, letter in enumerate(word) if letter == 'g']
count = len(positions)
print(f"Count: {count}, Positions: {positions}")