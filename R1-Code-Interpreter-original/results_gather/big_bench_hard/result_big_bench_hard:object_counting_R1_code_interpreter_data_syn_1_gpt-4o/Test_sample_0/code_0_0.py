# List of items
items = [
    "flute", "piano", "trombone", "stove", "stove", "stove", "stove",
    "violin", "accordion", "clarinet", "drum", "lamp", "lamp", "trumpet"
]

# List of musical instruments
musical_instruments = [
    "flute", "piano", "trombone", "violin", "accordion", "clarinet", "drum", "trumpet"
]

# Count the number of musical instruments
count = sum(item in musical_instruments for item in items)

print(count)