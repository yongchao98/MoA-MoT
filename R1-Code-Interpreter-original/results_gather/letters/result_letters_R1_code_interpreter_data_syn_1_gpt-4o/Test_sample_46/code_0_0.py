word = "hypophysectomies"
count = word.count('y')
positions = [i + 1 for i, letter in enumerate(word) if letter == 'y']
print(f"Count: {count}, Positions: {positions}")