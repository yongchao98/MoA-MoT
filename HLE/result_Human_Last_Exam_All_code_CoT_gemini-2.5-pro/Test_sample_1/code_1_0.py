import codecs

# Step 1: Define c1
# The reciprocal concept to logical "depth" is "crypticity".
# The third letter is 'y'.
c1 = "crypticity"[2]
print(f"c1 = '{c1}' (from 'crypticity')")

# Step 2: Define c2
# The missing word in Murray Gell-Man's quote is "options".
# The third character is 't'.
c2 = "options"[2]
print(f"c2 = '{c2}' (from 'options')")

# Step 3 & 4: Define c4
# The last author of the GELU paper is Gimpel, ending in 'l'.
# Rot13 of 'l' is 'y'.
c3 = "Gimpel"[-1]
c4 = codecs.encode(c3, 'rot_13')
print(f"c4 = '{c4}' (from Rot13 of '{c3}')")

# Step 5: Define c5
# Mars is closer in mass to the Moon.
# The second letter of "Moon" is 'o'.
c5 = "Moon"[1]
print(f"c5 = '{c5}' (from 'Moon')")

# Final step: Concatenate and make lowercase
result = (c1 + c2 + c4 + c5).lower()
print(f"\nFinal concatenation: {c1.lower()} + {c2.lower()} + {c4.lower()} + {c5.lower()} = {result}")
