penguins = [
    {"name": "Louis", "age": 7},
    {"name": "Bernard", "age": 5},
    {"name": "Vincent", "age": 9},
    {"name": "Gwen", "age": 8},
    {"name": "James", "age": 12}
]

oldest_penguin = max(penguins, key=lambda p: p["age"])
print(oldest_penguin["name"])