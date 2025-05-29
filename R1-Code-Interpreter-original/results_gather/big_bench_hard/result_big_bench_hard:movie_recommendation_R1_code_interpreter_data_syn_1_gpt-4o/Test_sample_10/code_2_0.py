# Define the movies and their attributes
movies = {
    "Dragonheart 3 The Sorcerer's Curse": {"genre": "Fantasy", "themes": ["Adventure", "Medieval"]},
    "Chi-Raq": {"genre": "Crime Drama", "themes": ["Social Issues", "Musical"]},
    "Stargate": {"genre": "Science Fiction", "themes": ["Adventure", "Ancient Civilizations"]},
    "Tae Guk Gi The Brotherhood of War": {"genre": "War Drama", "themes": ["Historical", "Brotherhood", "War"]},
}

# Define the target attributes based on the given movies
target_genres = ["Drama", "Historical"]
target_themes = ["Character Development", "Cultural Significance", "Historical"]

# Function to find the most similar movie
def find_similar_movie(movies, target_genres, target_themes):
    best_match = None
    max_matches = 0
    for movie, attributes in movies.items():
        genre_matches = sum(genre in attributes["genre"] for genre in target_genres)
        theme_matches = sum(theme in attributes["themes"] for theme in target_themes)
        total_matches = genre_matches + theme_matches
        if total_matches > max_matches:
            best_match = movie
            max_matches = total_matches
    return best_match

# Find and print the most similar movie
similar_movie = find_similar_movie(movies, target_genres, target_themes)
print(similar_movie)