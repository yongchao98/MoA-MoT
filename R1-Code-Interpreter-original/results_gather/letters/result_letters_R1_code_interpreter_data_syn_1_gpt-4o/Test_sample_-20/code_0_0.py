word = "pseudoparthenogenesis"
count = 0
positions = []

for index, char in enumerate(word):
    if char == 'g':
        count += 1
        positions.append(index + 1)

print(f"Count: {count}, Positions: {positions}")