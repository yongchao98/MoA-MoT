# Define a dictionary with movie options and their genres/themes
movies = {
    "Dead Man Walking": ["Drama", "Intense Storytelling"],
    "World's Greatest Dad": ["Dark Comedy", "Drama"],
    "The Bad and the Beautiful": ["Drama", "Film Industry"],
    "Hood of Horror": ["Horror", "Anthology"]
}

# Define the common themes/genres of the given movies
common_themes = ["Action", "Suspense", "Drama", "Intense Storytelling"]

# Function to find the most similar movie
def find_similar_movie(movies, common_themes):
    max_match = 0
    similar_movie = None
    for movie, themes in movies.items():
        match_count = len(set(themes) & set(common_themes))
        if match_count > max_match:
            max_match = match_count
            similar_movie = movie
    return similar_movie

# Find and print the most similar movie
similar_movie = find_similar_movie(movies, common_themes)
print(similar_movie)