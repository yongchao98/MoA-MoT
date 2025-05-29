# Define the genres of the given movies
given_movies_genres = {
    "Apollo 13": ["Drama", "History", "Adventure"],
    "Jurassic Park": ["Adventure", "Sci-Fi", "Thriller"],
    "Die Hard With a Vengeance": ["Action", "Thriller"],
    "Forrest Gump": ["Drama", "Romance", "Comedy"]
}

# Define the genres of the options
options_genres = {
    "Killer Movie": ["Horror", "Thriller"],
    "Stealth": ["Action", "Adventure", "Sci-Fi"],
    "The Last Man on Earth": ["Sci-Fi", "Horror"],
    "True Lies": ["Action", "Comedy", "Thriller"]
}

# Function to calculate similarity based on genre overlap
def calculate_similarity(given_genres, option_genres):
    similarity_score = 0
    for genre in option_genres:
        if genre in given_genres:
            similarity_score += 1
    return similarity_score

# Calculate similarity scores for each option
similarity_scores = {}
for option, genres in options_genres.items():
    total_similarity = 0
    for movie, movie_genres in given_movies_genres.items():
        total_similarity += calculate_similarity(movie_genres, genres)
    similarity_scores[option] = total_similarity

# Find the option with the highest similarity score
most_similar_movie = max(similarity_scores, key=similarity_scores.get)
print(most_similar_movie)