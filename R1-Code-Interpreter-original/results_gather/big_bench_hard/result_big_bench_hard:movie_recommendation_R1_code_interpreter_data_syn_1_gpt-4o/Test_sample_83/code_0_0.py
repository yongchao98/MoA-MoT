# Define the movies and their genres/themes
movies = {
    "Interview with the Vampire": ["gothic horror", "supernatural", "existentialism"],
    "The Shawshank Redemption": ["drama", "hope", "friendship", "redemption"],
    "Pulp Fiction": ["crime", "non-linear narrative", "dark humor", "stylized violence"],
    "The Usual Suspects": ["crime thriller", "deception", "twist ending"],
}

options = {
    "Caligula": ["historical drama", "controversial", "explicit"],
    "Anger Management": ["comedy", "therapy", "anger"],
    "The Fugitive": ["thriller", "suspense", "chase", "crime"],
    "Mommy": ["drama", "emotional", "mother-son relationship"],
}

# Function to find the most similar movie
def find_similar_movie(movies, options):
    # Count the number of matching genres/themes
    similarity_scores = {}
    for option, themes in options.items():
        score = 0
        for movie, movie_themes in movies.items():
            score += len(set(themes) & set(movie_themes))
        similarity_scores[option] = score
    # Find the option with the highest score
    most_similar = max(similarity_scores, key=similarity_scores.get)
    return most_similar

# Find and print the most similar movie
most_similar_movie = find_similar_movie(movies, options)
print(most_similar_movie)