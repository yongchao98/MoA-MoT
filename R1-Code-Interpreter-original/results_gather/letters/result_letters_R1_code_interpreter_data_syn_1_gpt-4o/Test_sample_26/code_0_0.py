word = "benzantialdoxime"
positions = [i + 1 for i, letter in enumerate(word) if letter == 'a']
count = len(positions)
print(f"Count: {count}, Positions: {positions}")