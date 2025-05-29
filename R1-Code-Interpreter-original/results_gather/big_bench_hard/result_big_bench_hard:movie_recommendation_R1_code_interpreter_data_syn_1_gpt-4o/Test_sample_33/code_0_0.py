# Define the characteristics of the given movies
given_movies = {
    "Braveheart": ["historical", "epic", "war", "heroism"],
    "Dances with Wolves": ["historical", "drama", "war", "transformation"],
    "Pulp Fiction": ["crime", "non-linear", "dark humor", "character development"],
    "Schindler's List": ["historical", "drama", "war", "heroism"]
}

# Define the characteristics of the options
options = {
    "Tickling Giants": ["documentary", "comedy", "political"],
    "Shadow of the Vampire": ["horror", "film-making", "vampire"],
    "The Usual Suspects": ["crime", "thriller", "complex narrative", "character development"],
    "On the Town": ["musical", "comedy", "romance"]
}

# Function to calculate similarity based on common characteristics
def calculate_similarity(movie_chars, option_chars):
    return len(set(movie_chars) & set(option_chars))

# Calculate similarity for each option
similarity_scores = {option: sum(calculate_similarity(chars, options[option]) for chars in given_movies.values()) for option in options}

# Find the option with the highest similarity score
most_similar_movie = max(similarity_scores, key=similarity_scores.get)

print(most_similar_movie)