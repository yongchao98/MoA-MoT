penguins = [
    {"name": "Louis", "age": 7},
    {"name": "Vincent", "age": 9},
    {"name": "Gwen", "age": 8}
]

count = sum(1 for penguin in penguins if penguin["age"] < 10)
print(count)