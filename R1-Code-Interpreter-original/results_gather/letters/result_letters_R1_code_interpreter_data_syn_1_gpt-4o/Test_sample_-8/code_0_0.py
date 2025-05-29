word = "counterclassifications"
positions = []
count = 0

for index, char in enumerate(word):
    if char == 's':
        positions.append(index + 1)  # Adding 1 to convert from 0-based to 1-based index
        count += 1

print(f"Count: {count}, Positions: {positions}")