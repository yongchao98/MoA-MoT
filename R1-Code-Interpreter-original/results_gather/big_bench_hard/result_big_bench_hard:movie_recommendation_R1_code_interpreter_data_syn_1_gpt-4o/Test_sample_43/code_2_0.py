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

# Define the most common genres to focus on
common_genres = ["Crime", "Drama", "Thriller"]

# Function to calculate similarity based on common genres
def calculate_similarity(option_genres):
    return len(set(common_genres) & set(option_genres))

# Calculate the total similarity for each option
similarity_scores = {option: calculate_similarity(option_genres) for option, option_genres in options.items()}

# Find the option with the highest similarity score
most_similar_option = max(similarity_scores, key=similarity_scores.get)
print(most_similar_option)