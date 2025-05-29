# Define the genres for each movie
movies = {
    "Jurassic Park": ["Science Fiction", "Adventure", "Thriller"],
    "The Silence of the Lambs": ["Psychological Horror", "Thriller"],
    "Schindler's List": ["Historical", "Drama"],
    "Braveheart": ["Historical", "War", "Drama"],
    "A Hard Day's Night": ["Musical", "Comedy"],
    "Showtime": ["Comedy"],
    "Forrest Gump": ["Drama", "Historical"],
    "Eddie the Eagle": ["Biographical", "Sports", "Comedy", "Drama"]
}

# Define the target genres based on the given movies
target_genres = ["Drama", "Historical", "Thriller"]

# Function to calculate similarity based on genre overlap
def calculate_similarity(movie_genres, target_genres):
    return len(set(movie_genres) & set(target_genres))

# Calculate similarity for each option
similarities = {movie: calculate_similarity(genres, target_genres) for movie, genres in movies.items()}

# Find the movie with the highest similarity
most_similar_movie = max(similarities, key=similarities.get)

print(most_similar_movie)