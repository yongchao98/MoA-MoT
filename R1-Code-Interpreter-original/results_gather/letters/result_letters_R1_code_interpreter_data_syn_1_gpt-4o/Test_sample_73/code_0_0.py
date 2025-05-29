word = 'zgmztnywyyovhlcdhindy'
count = 0
positions = []

for index, char in enumerate(word):
    if char == 'y':
        count += 1
        positions.append(index + 1)

print(f"Count: {count}, Positions: {positions}")