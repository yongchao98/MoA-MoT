# Define a dictionary with movie options and their genres/themes
movies = {
    "1984": ["Dystopian", "Science Fiction", "Drama"],
    "Good Night, and Good Luck": ["Historical", "Drama"],
    "Kiss of the Dragon": ["Action", "Martial Arts"],
    "District 9": ["Science Fiction", "Action", "Social Commentary"]
}

# Define the genres/themes of the original movies
original_movies = {
    "Pulp Fiction": ["Crime", "Drama", "Dark Humor"],
    "Forrest Gump": ["Drama", "Comedy", "Historical"],
    "Inside Out": ["Animation", "Drama", "Psychological"],
    "Edge of Tomorrow": ["Science Fiction", "Action", "Time Loop"]
}

# Function to find the most similar movie
def find_similar_movie(movies, original_movies):
    similarity_scores = {}
    for movie, genres in movies.items():
        score = 0
        for original, original_genres in original_movies.items():
            score += len(set(genres) & set(original_genres))
        similarity_scores[movie] = score
    # Find the movie with the highest similarity score
    similar_movie = max(similarity_scores, key=similarity_scores.get)
    return similar_movie

# Find and print the most similar movie
similar_movie = find_similar_movie(movies, original_movies)
print(similar_movie)