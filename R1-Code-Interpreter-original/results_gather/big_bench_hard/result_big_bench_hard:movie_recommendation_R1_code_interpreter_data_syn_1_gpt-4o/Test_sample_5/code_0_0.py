# Define the genres of the original films
original_genres = {
    "Heat": ["Crime", "Drama", "Thriller"],
    "The Fugitive": ["Action", "Thriller"],
    "Forrest Gump": ["Drama", "Comedy"],
    "The Silence of the Lambs": ["Thriller", "Crime", "Drama"]
}

# Define the genres of the options
options_genres = {
    "Death Race 2": ["Action", "Sci-Fi", "Thriller"],
    "Cannonball Run II": ["Comedy", "Action"],
    "Independence Day": ["Action", "Adventure", "Sci-Fi"],
    "Slumber Party Massacre II": ["Horror", "Comedy"]
}

# Function to calculate similarity based on genre overlap
def calculate_similarity(original_genres, option_genres):
    similarity_scores = {}
    for option, genres in option_genres.items():
        score = 0
        for original, original_genre_list in original_genres.items():
            overlap = set(genres).intersection(set(original_genre_list))
            score += len(overlap)
        similarity_scores[option] = score
    return similarity_scores

# Calculate similarity scores
similarity_scores = calculate_similarity(original_genres, options_genres)

# Find the option with the highest similarity score
most_similar_option = max(similarity_scores, key=similarity_scores.get)

print(most_similar_option)