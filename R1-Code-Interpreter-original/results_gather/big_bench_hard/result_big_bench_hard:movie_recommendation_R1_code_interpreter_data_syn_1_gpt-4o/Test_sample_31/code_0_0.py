# Define the genres for each movie
movies_genres = {
    "Batman": ["Action", "Crime", "Superhero"],
    "Pulp Fiction": ["Crime", "Drama", "Dark Comedy"],
    "Apollo 13": ["Historical", "Drama", "Adventure"],
    "Schindler's List": ["Historical", "Drama", "War"],
    "The Last Kiss": ["Romantic", "Drama"],
    "Europa": ["Science Fiction"],
    "Braveheart": ["Historical", "Drama", "War"],
    "Rise of the Planet of the Apes": ["Science Fiction", "Action"]
}

# Calculate similarity based on common genres
def calculate_similarity(target_movies, option):
    target_genres = set()
    for movie in target_movies:
        target_genres.update(movies_genres[movie])
    
    option_genres = set(movies_genres[option])
    common_genres = target_genres.intersection(option_genres)
    
    return len(common_genres)

# Target movies
target_movies = ["Batman", "Pulp Fiction", "Apollo 13", "Schindler's List"]

# Options
options = ["The Last Kiss", "Europa", "Braveheart", "Rise of the Planet of the Apes"]

# Find the most similar movie
most_similar = max(options, key=lambda option: calculate_similarity(target_movies, option))

print(most_similar)