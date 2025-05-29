# Define the family relationships
family = {
    "Benjamin": {"spouse": "Barbara", "children": ["Lucas"]},
    "Barbara": {"spouse": "Benjamin", "children": ["Lucas"]},
    "Lucas": {"spouse": "Olivia"},
    "Olivia": {"spouse": "Lucas"}
}

# Determine the relationship of Olivia to Lucas
relationship = family["Lucas"]["spouse"]

print(relationship)