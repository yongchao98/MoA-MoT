# Information about the artwork and its acquisition
artist = "Kurt Günther"
painting_title = "The Radionist"
creation_year = 1927
acquisition_year = 1967
acquiring_museum = "Galerie Neue Meister, part of the Staatliche Kunstsammlungen Dresden"

# Constructing the final answer string
output = (
    f"In {acquisition_year}, the museum that acquired the {creation_year} tempera painting "
    f"'{painting_title}' by {artist} was the {acquiring_museum}."
)

# Printing the result
print(output)