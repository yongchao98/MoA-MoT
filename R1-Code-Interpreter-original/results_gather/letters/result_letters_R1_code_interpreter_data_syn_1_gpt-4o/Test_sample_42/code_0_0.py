word = 'nonexpeditiousness'
positions = []

for index, char in enumerate(word):
    if char == 's':
        positions.append(index + 1)

count = len(positions)
print(f"Count: {count}, Positions: {positions}")