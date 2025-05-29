# Define the movies and their attributes
movies = {
    "Dragonheart 3 The Sorcerer's Curse": {"genre": "Fantasy", "themes": ["Adventure", "Medieval"]},
    "Chi-Raq": {"genre": "Crime Drama", "themes": ["Social Issues", "Musical"]},
    "Stargate": {"genre": "Science Fiction", "themes": ["Adventure", "Ancient Civilizations"]},
    "Tae Guk Gi The Brotherhood of War": {"genre": "War Drama", "themes": ["Historical", "Brotherhood", "War"]},
}

# Define the target attributes based on the given movies
target_genres = ["Drama", "Historical"]
target_themes = ["Character Development", "Cultural Significance"]

# Function to find the most similar movie
def find_similar_movie(movies, target_genres, target_themes):
    for movie, attributes in movies.items():
        if any(genre in attributes["genre"] for genre in target_genres) and any(theme in attributes["themes"] for theme in target_themes):
            return movie
    return None

# Find and print the most similar movie
similar_movie = find_similar_movie(movies, target_genres, target_themes)
print(similar_movie)