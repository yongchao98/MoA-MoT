# Define the genres of the options
options_genres = {
    "Pleasantville": ["Fantasy", "Comedy", "Drama"],
    "The Mask": ["Superhero", "Comedy"],
    "To Catch a Thief": ["Romantic", "Thriller"],
    "Microwave Massacre": ["Comedy", "Horror"]
}

# Define the common genres of the given movies
common_genres = ["Action", "Adventure", "Drama", "Science Fiction", "Superhero"]

# Find the option with the most similar genres
similarity_scores = {}
for option, genres in options_genres.items():
    similarity_scores[option] = len(set(genres) & set(common_genres))

# Find the option with the highest similarity score
most_similar_movie = max(similarity_scores, key=similarity_scores.get)
print(most_similar_movie)