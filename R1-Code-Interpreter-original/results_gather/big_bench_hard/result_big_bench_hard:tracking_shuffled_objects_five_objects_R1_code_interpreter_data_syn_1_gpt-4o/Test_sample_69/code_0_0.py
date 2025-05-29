# Initial pairings
partners = {
    "Alice": "Karl",
    "Bob": "Melissa",
    "Claire": "Lola",
    "Dave": "Sam",
    "Eve": "Izzi"
}

# Apply the partner switches
# 1. Claire and Dave switch partners
partners["Claire"], partners["Dave"] = partners["Dave"], partners["Claire"]

# 2. Eve and Bob switch partners
partners["Eve"], partners["Bob"] = partners["Bob"], partners["Eve"]

# 3. Claire and Alice switch partners
partners["Claire"], partners["Alice"] = partners["Alice"], partners["Claire"]

# 4. Eve and Alice switch partners
partners["Eve"], partners["Alice"] = partners["Alice"], partners["Eve"]

# 5. Claire and Dave switch partners
partners["Claire"], partners["Dave"] = partners["Dave"], partners["Claire"]

# Find out who Alice is dancing with
print(partners["Alice"])