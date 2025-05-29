# Define the movies and their genres/themes
movies = {
    "Braveheart": ["historical", "epic", "war", "heroism"],
    "Dances with Wolves": ["historical", "drama", "cultural", "war"],
    "The Shawshank Redemption": ["drama", "hope", "friendship", "redemption"],
    "The Usual Suspects": ["crime", "thriller", "mystery", "deception"],
    "3000 Miles to Graceland": ["crime", "action", "heist", "betrayal"],
    "Crimson Tide": ["thriller", "drama", "military", "leadership"],
    "Best Men": ["comedy", "crime", "friendship", "heist"],
    "A Very Harold & Kumar 3D Christmas": ["comedy", "friendship", "holiday", "adventure"]
}

# Calculate similarity based on common themes
def calculate_similarity(target_movie, options):
    target_themes = set(movies[target_movie])
    similarities = {}
    for option in options:
        option_themes = set(movies[option])
        common_themes = target_themes.intersection(option_themes)
        similarities[option] = len(common_themes)
    return similarities

# Define the target movies and options
target_movies = ["Braveheart", "Dances with Wolves", "The Shawshank Redemption", "The Usual Suspects"]
options = ["3000 Miles to Graceland", "Crimson Tide", "Best Men", "A Very Harold & Kumar 3D Christmas"]

# Calculate total similarity for each option
total_similarity = {option: 0 for option in options}
for target in target_movies:
    similarity = calculate_similarity(target, options)
    for option, score in similarity.items():
        total_similarity[option] += score

# Find the option with the highest similarity score
most_similar_movie = max(total_similarity, key=total_similarity.get)
print(most_similar_movie)