# Define the genres of the given movies
movies_genres = {
    "Jurassic Park": ["Science Fiction", "Adventure", "Thriller", "Action"],
    "The Silence of the Lambs": ["Psychological Horror", "Thriller"],
    "Batman": ["Superhero", "Action", "Thriller"],
    "American Beauty": ["Drama", "Dark Comedy", "Psychological"],
}

# Define the genres of the options
options_genres = {
    "Night of the Living Dead": ["Horror", "Thriller"],
    "Forrest Gump": ["Drama", "Romance", "Comedy"],
    "Ghost Ship": ["Horror", "Thriller", "Psychological"],
    "Revenge for Jolly!": ["Comedy", "Drama"],
}

# Function to calculate similarity based on common genres
def calculate_similarity(movie_genres, option_genres):
    return len(set(movie_genres) & set(option_genres))

# Calculate similarity for each option
similarity_scores = {option: 0 for option in options_genres}
for movie, genres in movies_genres.items():
    for option, option_genres in options_genres.items():
        similarity_scores[option] += calculate_similarity(genres, option_genres)

# Find the option with the highest similarity score
most_similar_movie = max(similarity_scores, key=similarity_scores.get)
print(most_similar_movie)