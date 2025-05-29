# Define the genres and themes of the given movies
movies = {
    "Léon: The Professional": ["crime", "thriller", "action"],
    "The Matrix": ["science fiction", "action"],
    "Forrest Gump": ["drama", "romance", "comedy"],
    "Pulp Fiction": ["crime", "drama", "non-linear narrative"]
}

# Define the genres and themes of the options
options = {
    "Romper Stomper": ["drama", "crime", "social issues"],
    "Fame": ["musical", "drama"],
    "Up": ["animated", "adventure", "friendship"],
    "The Browning Version": ["drama"]
}

# Function to find the most similar movie
def find_similar_movie(movies, options):
    # Calculate similarity based on common genres/themes
    similarity_scores = {}
    for option, option_genres in options.items():
        score = 0
        for movie_genres in movies.values():
            score += len(set(option_genres) & set(movie_genres))
        similarity_scores[option] = score
    # Find the option with the highest similarity score
    most_similar = max(similarity_scores, key=similarity_scores.get)
    return most_similar

# Find and print the most similar movie
most_similar_movie = find_similar_movie(movies, options)
print(most_similar_movie)