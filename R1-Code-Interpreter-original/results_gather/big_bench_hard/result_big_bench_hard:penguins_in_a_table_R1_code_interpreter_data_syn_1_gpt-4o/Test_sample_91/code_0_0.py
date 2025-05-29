# List of penguins with their ages
penguins = [
    {"name": "Louis", "age": 7},
    {"name": "Bernard", "age": 5},
    {"name": "Vincent", "age": 9},
    {"name": "Gwen", "age": 8},
    {"name": "James", "age": 12}
]

# Count penguins with age less than 10
count_under_10 = sum(1 for penguin in penguins if penguin["age"] < 10)

# Output the result
print(count_under_10)