# Initial pairings
partners = {
    "Alice": "Melissa",
    "Bob": "Lola",
    "Claire": "Jamie",
    "Dave": "Patrick",
    "Eve": "Ophelia"
}

# Step 1: Bob and Dave switch partners
partners["Bob"], partners["Dave"] = partners["Dave"], partners["Bob"]

# Step 2: Claire and Eve switch partners
partners["Claire"], partners["Eve"] = partners["Eve"], partners["Claire"]

# Step 3: Dave and Alice switch partners
partners["Dave"], partners["Alice"] = partners["Alice"], partners["Dave"]

# Step 4: Claire and Eve switch partners again
partners["Claire"], partners["Eve"] = partners["Eve"], partners["Claire"]

# Step 5: Eve and Bob switch partners
partners["Eve"], partners["Bob"] = partners["Bob"], partners["Eve"]

# Output the final partner of Bob
print(partners["Bob"])