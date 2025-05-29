word = 'mbxxoatfbjtfpndbkyssmm'
count = 0
positions = []

for index, char in enumerate(word):
    if char == 'm':
        count += 1
        positions.append(index + 1)  # Convert 0-based index to 1-based

print(f"Count: {count}, Positions: {positions}")