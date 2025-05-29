items = ["piano", "trombone", "clarinet", "goat", "accordion", "trumpet"]
musical_instruments = ["piano", "trombone", "clarinet", "accordion", "trumpet"]

count = sum(1 for item in items if item in musical_instruments)
print(count)