word = "coumarilic"
count = 0
positions = []

for index, char in enumerate(word):
    if char == 'm':
        count += 1
        positions.append(index + 1)

print(f"Count: {count}, Positions: {positions}")