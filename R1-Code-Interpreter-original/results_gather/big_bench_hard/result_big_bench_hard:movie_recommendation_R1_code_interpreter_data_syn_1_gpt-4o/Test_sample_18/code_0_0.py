# Define the characteristics of the given movies
given_movies = {
    "Pulp Fiction": ["crime", "non-linear", "dark humor", "neo-noir"],
    "The Shawshank Redemption": ["drama", "hope", "friendship", "prison"],
    "Schindler's List": ["historical", "drama", "Holocaust", "emotional"],
    "Dances with Wolves": ["historical", "epic", "cultural", "landscapes"]
}

# Define the characteristics of the options
options = {
    "Easter Parade": ["musical", "dance", "romance"],
    "Captain America The First Avenger": ["superhero", "action", "adventure"],
    "Braveheart": ["historical", "epic", "drama", "emotional"],
    "Manic": ["drama", "teen", "institution"]
}

# Function to calculate similarity based on common characteristics
def calculate_similarity(movie_chars, option_chars):
    return len(set(movie_chars) & set(option_chars))

# Calculate similarity for each option
similarity_scores = {option: sum(calculate_similarity(chars, options[option]) for chars in given_movies.values()) for option in options}

# Find the option with the highest similarity score
most_similar_movie = max(similarity_scores, key=similarity_scores.get)

print(most_similar_movie)