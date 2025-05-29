# Define a dictionary with movie options and their genres/themes
movies = {
    "Frankenstein Unbound": ["Science Fiction", "Horror", "Time Travel"],
    "Jon Stewart Has Left the Building": ["Comedy", "Documentary"],
    "Clear and Present Danger": ["Action", "Thriller", "Drama", "Espionage"],
    "Little Giants": ["Family", "Sports", "Comedy"]
}

# Define the target themes based on the given movies
target_themes = ["Drama", "Thriller", "Action", "Survival", "Redemption"]

# Function to find the most similar movie based on matching themes
def find_similar_movie(movies, target_themes):
    max_matches = 0
    similar_movie = None
    for movie, themes in movies.items():
        matches = len(set(themes) & set(target_themes))
        if matches > max_matches:
            max_matches = matches
            similar_movie = movie
    return similar_movie

# Find and print the most similar movie
similar_movie = find_similar_movie(movies, target_themes)
print(similar_movie)