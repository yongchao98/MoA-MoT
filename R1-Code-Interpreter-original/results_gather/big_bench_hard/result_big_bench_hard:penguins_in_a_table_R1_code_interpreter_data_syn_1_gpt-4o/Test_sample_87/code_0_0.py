# List of penguins with their ages after removing Bernard
penguins = [
    {"name": "Louis", "age": 7},
    {"name": "Vincent", "age": 9},
    {"name": "Gwen", "age": 8}
]

# Count penguins older than 5 years
count = sum(1 for penguin in penguins if penguin["age"] > 5)

print(count)