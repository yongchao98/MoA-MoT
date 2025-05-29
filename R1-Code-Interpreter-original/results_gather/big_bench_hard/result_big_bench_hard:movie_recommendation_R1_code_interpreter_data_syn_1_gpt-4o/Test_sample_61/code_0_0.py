# Define a dictionary with movie options and their genres/themes
movies = {
    "Boy": ["coming-of-age", "drama"],
    "Swiss Family Robinson": ["adventure", "family"],
    "And Then There Were None": ["mystery", "thriller", "crime"],
    "Crimson Tide": ["thriller", "drama", "military", "suspense"]
}

# Define the genres/themes of the original movies
original_movies_themes = ["drama", "thriller", "suspense", "freedom", "justice"]

# Function to find the most similar movie
def find_similar_movie(movies, original_themes):
    max_overlap = 0
    best_match = None
    for movie, themes in movies.items():
        overlap = len(set(themes) & set(original_themes))
        if overlap > max_overlap:
            max_overlap = overlap
            best_match = movie
    return best_match

# Find and print the most similar movie
similar_movie = find_similar_movie(movies, original_movies_themes)
print(similar_movie)