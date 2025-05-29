# Define the movies and their genres/themes
movies = {
    "The Sixth Sense": ["psychological thriller", "supernatural", "twist ending"],
    "The Matrix": ["science fiction", "action", "philosophical"],
    "Forrest Gump": ["drama", "romance", "comedy", "historical"],
    "The Shawshank Redemption": ["drama", "hope", "friendship"],
    "Street Fighter II The Animated Movie": ["animated", "action", "martial arts"],
    "The Sheltering Sky": ["drama", "existential", "human condition"],
    "The Boy Who Could Fly": ["fantasy", "drama", "imagination"],
    "Terminator 2 Judgment Day": ["science fiction", "action", "artificial intelligence", "emotional"]
}

# Calculate similarity based on common genres/themes
def calculate_similarity(target_movies, option_movies):
    similarity_scores = {}
    for option, themes in option_movies.items():
        similarity_scores[option] = sum(1 for theme in themes if any(theme in target for target in target_movies.values()))
    return similarity_scores

# Target movies
target_movies = {key: movies[key] for key in ["The Sixth Sense", "The Matrix", "Forrest Gump", "The Shawshank Redemption"]}

# Option movies
option_movies = {key: movies[key] for key in ["Street Fighter II The Animated Movie", "The Sheltering Sky", "The Boy Who Could Fly", "Terminator 2 Judgment Day"]}

# Calculate and print similarity scores
similarity_scores = calculate_similarity(target_movies, option_movies)
print(similarity_scores)