word = 'rmzfrwkdjocilzylmvmgayms'
count = 0
positions = []

for index, char in enumerate(word):
    if char == 'm':
        count += 1
        positions.append(index + 1)  # Adding 1 to convert 0-based index to 1-based position

print(f"Count: {count}, Positions: {positions}")