# Initial pairings
pairings = {
    "Alice": "Ophelia",
    "Bob": "Jamie",
    "Claire": "Sam",
    "Dave": "Melissa",
    "Eve": "Rodrigo"
}

# Apply the partner switches
# Eve and Dave switch partners
pairings["Dave"], pairings["Eve"] = pairings["Eve"], pairings["Dave"]

# Alice and Bob switch partners
pairings["Alice"], pairings["Bob"] = pairings["Bob"], pairings["Alice"]

# Dave and Bob switch partners
pairings["Dave"], pairings["Bob"] = pairings["Bob"], pairings["Dave"]

# Bob and Alice switch partners
pairings["Bob"], pairings["Alice"] = pairings["Alice"], pairings["Bob"]

# Alice and Claire switch partners
pairings["Alice"], pairings["Claire"] = pairings["Claire"], pairings["Alice"]

# Determine who Bob is dancing with at the end
final_partner_bob = pairings["Bob"]
print(final_partner_bob)