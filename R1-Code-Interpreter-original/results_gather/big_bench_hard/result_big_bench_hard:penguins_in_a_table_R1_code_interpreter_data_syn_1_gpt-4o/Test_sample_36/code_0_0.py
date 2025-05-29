penguins = {
    "Louis": 7,
    "Vincent": 9,
    "Gwen": 8,
    "James": 12
}

oldest_penguin = max(penguins, key=penguins.get)
print(oldest_penguin)