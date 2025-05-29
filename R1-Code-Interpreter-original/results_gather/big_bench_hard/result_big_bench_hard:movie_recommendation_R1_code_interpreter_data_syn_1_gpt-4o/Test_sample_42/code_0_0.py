# Define the themes of the given movies
given_movies_themes = {
    "Dances with Wolves": ["historical", "drama", "cultural exchange", "transformation"],
    "The Lion King": ["animated", "leadership", "redemption", "circle of life"],
    "The Fugitive": ["thriller", "justice", "pursuit", "survival"],
    "The Shawshank Redemption": ["drama", "hope", "friendship", "redemption"]
}

# Define the themes of the options
options_themes = {
    "Schindler's List": ["historical", "drama", "redemption", "humanity"],
    "Day for Night": ["filmmaking", "relationships", "drama"],
    "Man of the House": ["comedy", "protection", "humor"],
    "The Neon Demon": ["psychological", "horror", "beauty", "obsession"]
}

# Function to calculate similarity based on common themes
def calculate_similarity(given_themes, option_themes):
    return len(set(given_themes) & set(option_themes))

# Calculate similarity for each option
similarity_scores = {option: sum(calculate_similarity(themes, options_themes[option]) for themes in given_movies_themes.values()) for option in options_themes}

# Find the option with the highest similarity score
most_similar_movie = max(similarity_scores, key=similarity_scores.get)

print(most_similar_movie)