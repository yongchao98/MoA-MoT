penguins = {
    "Louis": 7,
    "Bernard": 5,
    "Vincent": 9,
    "Gwen": 8
}

oldest_penguin = max(penguins, key=penguins.get)
print(oldest_penguin)