# The sound elements to be ranked
sound_elements = {
    'Music': 1,
    'Sound Effects': 2,
    'Speech': 3
}

# Sort the elements based on their rank (value)
ranked_elements = sorted(sound_elements.keys(), key=lambda item: sound_elements[item])

# Join the sorted list into a comma-separated string for display
final_ranking = ", ".join(ranked_elements)

# Print the final hierarchical order
print(final_ranking)