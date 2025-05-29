# Define the genres and themes of the given movies
given_movies = {
    "Pulp Fiction": ["crime", "dark humor", "non-linear"],
    "Forrest Gump": ["drama", "comedy", "historical"],
    "The Fugitive": ["action", "thriller", "innocence"],
    "Braveheart": ["historical", "epic", "freedom"]
}

# Define the genres and themes of the options
options = {
    "Little Women": ["drama", "period", "family"],
    "Dances with Wolves": ["historical", "epic", "cultural"],
    "Bring It On All or Nothing": ["comedy", "teen", "sports"],
    "Dead Alive": ["horror", "comedy", "gore"]
}

# Function to calculate similarity based on common themes
def calculate_similarity(movie_themes, option_themes):
    return len(set(movie_themes) & set(option_themes))

# Calculate similarity for each option
similarity_scores = {option: sum(calculate_similarity(themes, options[option]) for themes in given_movies.values()) for option in options}

# Find the option with the highest similarity score
most_similar_movie = max(similarity_scores, key=similarity_scores.get)

print(most_similar_movie)