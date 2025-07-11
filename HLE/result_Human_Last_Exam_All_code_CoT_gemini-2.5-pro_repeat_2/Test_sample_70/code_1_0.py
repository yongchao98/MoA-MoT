import collections

def get_flag_rank(description):
    """
    Determines the algebraic rank of a flag based on its structural description.
    This is a simplified model based on geometric patterns.

    - Rank 1: Simple horizontal or vertical stripes without other features.
    - Rank 2: A base of stripes (rank 1) with a simple, centered emblem,
              or intersecting uniform stripes.
    - Rank > 2: More complex patterns like multiple emblems, complex shapes (e.g., pall),
                or cantons with detailed elements.
    """
    if "horizontal stripes with central emblem" in description:
        return 2
    if "vertical stripes with central emblem" in description:
        return 2
    if "horizontal stripes intersected by vertical stripe" in description:
        return 2
    if "cross" in description and "complex" not in description:
        return 2
    if "horizontal stripes" in description or "vertical stripes" in description:
        return 1
    # Default for more complex or unclassified flags
    return 3

def find_flags_with_same_rank_as_denmark():
    """
    Finds African flags with the same algebraic rank as the flag of Denmark.
    """
    # 1. Determine the rank of the Danish flag
    denmark_description = "red field with a white cross"
    denmark_rank = get_flag_rank(denmark_description)
    print(f"The analysis begins with the flag of Denmark.")
    print(f"Its structure of a simple cross gives it a maximal linear algebraic rank of {denmark_rank}.")
    print("-" * 30)
    print(f"Searching for African flags with a rank of {denmark_rank}...\n")

    # A dictionary of African nations and their flag descriptions for our model
    african_flags = {
        "Algeria": "vertical stripes with central emblem (crescent and star)",
        "Angola": "horizontal stripes with central emblem (machete and gear)",
        "Benin": "vertical stripes next to horizontal stripes", # Rank > 2
        "Botswana": "horizontal stripes",
        "Burundi": "complex cross (saltire) with central emblem", # Rank > 2
        "Cape Verde": "horizontal stripes with complex emblem (circle of stars, off-center)", # Rank > 2
        "Central African Republic": "horizontal stripes intersected by vertical stripe",
        "Egypt": "horizontal stripes with central emblem (Eagle of Saladin)",
        "Eswatini": "horizontal stripes with central emblem (shield and spears)",
        "Ethiopia": "horizontal stripes with central emblem (star pentagram)",
        "Ghana": "horizontal stripes with central emblem (black star)",
        "Kenya": "horizontal stripes with central emblem (Maasai shield and spears)",
        "Libya": "horizontal stripes with central emblem (crescent and star)",
        "Mali": "vertical stripes",
        "Mozambique": "horizontal stripes with complex emblem on hoist triangle", # Rank > 2
        "Niger": "horizontal stripes with central emblem (orange disc)",
        "Nigeria": "vertical stripes",
        "Somalia": "solid field with central emblem (star)", # Rank > 2 due to star shape
        "South Africa": "complex pall design", # Rank > 2
        "Tunisia": "solid field with complex central emblem (disc with crescent and star)", # Rank > 2
        "Uganda": "horizontal stripes with central emblem (grey crowned crane)", # Crane shape is complex, rank > 2
    }

    matching_countries = []
    for country, desc in sorted(african_flags.items()):
        rank = get_flag_rank(desc)
        if rank == denmark_rank:
            matching_countries.append(country)

    print("The following African nations have flags with the same rank as Denmark:")
    for country in matching_countries:
        print(f"- {country}")

# Execute the function
find_flags_with_same_rank_as_denmark()

# The final equation demonstrates the principle for a flag with horizontal stripes and a central emblem.
# Let v_ones be a vector of 1s. Let v_emblem be a vector representing the emblem's shape.
# Any row in the flag matrix (R_i) can be expressed as: R_i = c_i * v_ones + k * v_emblem
# where c_i is the color of the stripe and k is a constant (0 for rows without the emblem).
# The row space is spanned by {v_ones, v_emblem}, so the rank is 2.
print("\nFinal equation principle: Rank = dimension(span({vector_for_stripes, vector_for_emblem})) = 2")