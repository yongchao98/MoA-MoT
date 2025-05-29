# Initial partners
partners = {
    "Alice": "Rodrigo",
    "Bob": "Sam",
    "Claire": "Lola"
}

# First switch: Claire and Bob switch partners
partners["Claire"], partners["Bob"] = partners["Bob"], partners["Claire"]

# Second switch: Bob and Alice switch partners
partners["Bob"], partners["Alice"] = partners["Alice"], partners["Bob"]

# Third switch: Claire and Alice switch partners
partners["Claire"], partners["Alice"] = partners["Alice"], partners["Claire"]

# Final partner of Claire
print(partners["Claire"])