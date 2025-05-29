# Define the data for penguins and giraffes
penguins = [
    {"name": "Louis", "age": 7},
    {"name": "Bernard", "age": 5},
    {"name": "Vincent", "age": 9},
    {"name": "Gwen", "age": 8}
]

giraffes = [
    {"name": "Jody", "age": 5},
    {"name": "Gladys", "age": 10},
    {"name": "Marian", "age": 2},
    {"name": "Donna", "age": 9}
]

# Count animals more than 5 years old
count = 0

# Check penguins
for penguin in penguins:
    if penguin["age"] > 5:
        count += 1

# Check giraffes
for giraffe in giraffes:
    if giraffe["age"] > 5:
        count += 1

# Output the result
print(count)