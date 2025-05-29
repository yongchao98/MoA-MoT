# Define the genres of the given movies
movies = {
    "Schindler's List": ["Historical", "Drama"],
    "Jurassic Park": ["Science Fiction", "Adventure"],
    "The Silence of the Lambs": ["Psychological", "Horror", "Thriller"],
    "Forrest Gump": ["Drama", "Historical"]
}

# Define the genres of the options
options = {
    "Batman": ["Action", "Superhero"],
    "Alien Resurrection": ["Science Fiction", "Horror"],
    "A Tale of Two Cities": ["Historical", "Drama"],
    "The Quiet American": ["Drama", "Historical"]
}

# Function to find the most similar movie
def find_similar_movie(movies, options):
    # Count the number of matching genres
    similarity_scores = {}
    for option, option_genres in options.items():
        score = 0
        for movie_genres in movies.values():
            score += len(set(option_genres) & set(movie_genres))
        similarity_scores[option] = score
    # Find the option with the highest score
    most_similar = max(similarity_scores, key=similarity_scores.get)
    return most_similar

# Find and print the most similar movie
most_similar_movie = find_similar_movie(movies, options)
print(most_similar_movie)