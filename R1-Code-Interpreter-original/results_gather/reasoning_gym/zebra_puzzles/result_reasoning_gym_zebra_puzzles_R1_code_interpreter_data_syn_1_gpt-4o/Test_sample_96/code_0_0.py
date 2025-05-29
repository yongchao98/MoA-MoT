# Assigning the known values to the houses
houses = {
    1: {"name": "Bob", "cigar": "Prince"},
    2: {"name": "Alice", "phone": "Google Pixel 6", "animal": "bird"},
    3: {"name": "Carol", "cigar": "Dunhill", "animal": "fish", "phone": "OnePlus 9"},
    4: {"name": "Arnold", "animal": "dog", "phone": "Huawei P50", "cigar": "Pall Mall"}
}

# Extracting the name of the person in House 1
name_house_1 = houses[1]["name"]
print(name_house_1)