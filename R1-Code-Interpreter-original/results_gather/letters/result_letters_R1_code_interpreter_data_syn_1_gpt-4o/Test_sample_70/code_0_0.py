word = "tychoparthenogenesis"
count = 0
positions = []

for index, char in enumerate(word):
    if char == 'y':
        count += 1
        positions.append(index + 1)  # Adding 1 to convert from 0-based to 1-based index

print(f"Count: {count}, Positions: {positions}")