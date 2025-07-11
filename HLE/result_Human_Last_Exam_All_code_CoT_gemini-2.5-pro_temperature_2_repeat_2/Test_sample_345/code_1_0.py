# The task is to identify the museum that acquired the 1927 painting "The Radionist" by Kurt Günther in 1967.
# Research indicates that the Staatliche Kunstsammlungen Dresden (Dresden State Art Collections) acquired the painting in that year.

# Define variables for the key information.
artist = "Kurt Günther"
painting_title = "The Radionist"
creation_year = 1927
acquisition_year = 1967
museum = "Staatliche Kunstsammlungen Dresden"

# Create a formatted string to present the answer, including the required numbers.
# This answers the question of which museum acquired the painting in the specified year.
answer_text = (
    f"The museum that acquired the {creation_year} tempera painting '{painting_title}' "
    f"by {artist} in {acquisition_year} was the {museum}."
)

# Print the final answer sentence.
print(answer_text)