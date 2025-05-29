# Initial pairings
partners = {
    "Alice": "Sam",
    "Bob": "Izzi",
    "Claire": "Jamie",
    "Dave": "Ophelia",
    "Eve": "Patrick"
}

# Step 1: Bob and Dave switch partners
partners["Bob"], partners["Dave"] = partners["Dave"], partners["Bob"]

# Step 2: Alice and Eve switch partners
partners["Alice"], partners["Eve"] = partners["Eve"], partners["Alice"]

# Step 3: Claire and Dave switch partners
partners["Claire"], partners["Dave"] = partners["Dave"], partners["Claire"]

# Step 4: Bob and Eve switch partners
partners["Bob"], partners["Eve"] = partners["Eve"], partners["Bob"]

# Step 5: Alice and Bob switch partners
partners["Alice"], partners["Bob"] = partners["Bob"], partners["Alice"]

# Final partner of Bob
print(partners["Bob"])