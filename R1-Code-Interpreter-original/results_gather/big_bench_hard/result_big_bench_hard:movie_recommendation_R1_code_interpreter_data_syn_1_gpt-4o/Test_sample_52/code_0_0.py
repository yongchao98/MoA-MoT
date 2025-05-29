# Define the genres and themes of the given movies and options
movies = {
    "Braveheart": ["historical", "epic", "war", "freedom", "rebellion", "leadership"],
    "The Usual Suspects": ["crime", "thriller", "mystery", "suspense"],
    "Crimson Tide": ["military", "thriller", "leadership", "conflict", "tension"],
    "The Fugitive": ["action", "thriller", "pursuit", "justice", "survival"],
}

options = {
    "Get Shorty": ["crime", "comedy", "mobsters", "film industry"],
    "The Firm": ["legal", "thriller", "conspiracy", "tension", "suspense"],
    "Pillow Talk": ["romantic", "comedy", "light-hearted"],
    "Hairspray": ["musical", "comedy-drama", "social change", "acceptance"],
}

# Function to find the most similar movie
def find_similar_movie(movies, options):
    # Calculate similarity based on common themes and genres
    similarity_scores = {}
    for option, option_genres in options.items():
        score = 0
        for movie, movie_genres in movies.items():
            common_elements = set(option_genres) & set(movie_genres)
            score += len(common_elements)
        similarity_scores[option] = score
    # Find the option with the highest similarity score
    most_similar = max(similarity_scores, key=similarity_scores.get)
    return most_similar

# Find and print the most similar movie
similar_movie = find_similar_movie(movies, options)
print(similar_movie)