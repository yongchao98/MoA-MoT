# Define the genres for each movie
given_movies = {
    "The Usual Suspects": ["Crime", "Thriller"],
    "Interview with the Vampire": ["Horror", "Drama"],
    "The Shawshank Redemption": ["Drama"],
    "Pulp Fiction": ["Crime", "Drama"]
}

options = {
    "Toy Soldiers": ["Action", "Drama"],
    "The Fugitive": ["Crime", "Thriller"],
    "The Wasp Woman": ["Horror", "Sci-Fi"],
    "Baxter": ["Horror", "Drama"]
}

# Function to calculate similarity based on genre
def calculate_similarity(movie_genres, option_genres):
    return len(set(movie_genres) & set(option_genres))

# Calculate the total similarity for each option
similarity_scores = {option: 0 for option in options}

for movie, genres in given_movies.items():
    for option, option_genres in options.items():
        similarity_scores[option] += calculate_similarity(genres, option_genres)

# Find the option with the highest similarity score
most_similar_option = max(similarity_scores, key=similarity_scores.get)
print(most_similar_option)