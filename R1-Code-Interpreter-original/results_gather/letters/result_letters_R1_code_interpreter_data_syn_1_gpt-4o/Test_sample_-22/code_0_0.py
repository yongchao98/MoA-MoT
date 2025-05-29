word = "palaeoclimatological"
count = 0
positions = []

for index, letter in enumerate(word):
    if letter == 'a':
        count += 1
        positions.append(index + 1)  # Adding 1 to convert from 0-based to 1-based index

print(f"Count: {count}, Positions: {positions}")