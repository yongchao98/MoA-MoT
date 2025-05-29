penguins = [
    {"name": "Louis", "age": 7, "weight": 11},
    {"name": "Bernard", "age": 5, "weight": 13},
    {"name": "Vincent", "age": 9, "weight": 11},
    {"name": "Gwen", "age": 8, "weight": 15},
    {"name": "James", "age": 12, "weight": 12}
]

count = sum(1 for penguin in penguins if penguin["age"] < 8 and penguin["weight"] < 12)
print(count)