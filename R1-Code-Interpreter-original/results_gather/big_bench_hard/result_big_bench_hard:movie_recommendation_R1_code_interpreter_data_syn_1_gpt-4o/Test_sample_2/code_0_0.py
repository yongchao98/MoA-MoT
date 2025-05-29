# Define a dictionary with movie options and their genres/themes
movies = {
    "Orlando": ["Historical", "Drama", "Fantasy"],
    "Guilty of Romance": ["Crime", "Drama", "Erotic"],
    "Forrest Gump": ["Drama", "Comedy", "Historical"],
    "All the Real Girls": ["Romance", "Drama", "Indie"]
}

# Define the common themes/genres of the given movies
common_themes = ["Historical", "Drama", "Character-driven"]

# Function to find the most similar movie
def find_similar_movie(movies, common_themes):
    similarity_scores = {}
    for movie, themes in movies.items():
        # Calculate similarity score based on common themes
        similarity_scores[movie] = len(set(themes) & set(common_themes))
    
    # Find the movie with the highest similarity score
    similar_movie = max(similarity_scores, key=similarity_scores.get)
    return similar_movie

# Find and print the most similar movie
similar_movie = find_similar_movie(movies, common_themes)
print(similar_movie)