# Define a dictionary with movie titles and their genres
movies_genres = {
    "Forrest Gump": ["Drama", "Romance", "Comedy"],
    "The Fugitive": ["Action", "Thriller"],
    "Independence Day": ["Science Fiction", "Action"],
    "Fargo": ["Crime", "Drama", "Dark Comedy"],
    "Comic-Con Episode IV A Fan's Hope": ["Documentary"],
    "Harley Davidson and the Marlboro Man": ["Action", "Crime", "Drama"],
    "Heat": ["Crime", "Drama", "Thriller"],
    "Scratch": ["Documentary"]
}

# Define the target genres based on the given movies
target_genres = ["Drama", "Action", "Crime", "Thriller"]

# Function to find the most similar movie
def find_similar_movie(options, target_genres):
    similarity_scores = {}
    for movie in options:
        genres = movies_genres[movie]
        # Calculate similarity score based on common genres
        similarity_score = len(set(genres) & set(target_genres))
        similarity_scores[movie] = similarity_score
    # Find the movie with the highest similarity score
    most_similar_movie = max(similarity_scores, key=similarity_scores.get)
    return most_similar_movie

# Options to consider
options = ["Comic-Con Episode IV A Fan's Hope", "Harley Davidson and the Marlboro Man", "Heat", "Scratch"]

# Find and print the most similar movie
most_similar_movie = find_similar_movie(options, target_genres)
print(most_similar_movie)