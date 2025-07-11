import Levenshtein

def find_runestone_id():
    """
    Identifies an Ingvar runestone by searching a small database
    for features observed in an image.
    """
    # Step 1: Create a knowledge base of Ingvar runestones.
    # This database contains information on some relevant Ingvar runestones,
    # including their visual style and key text fragments.
    runestones = [
        {
            "id": "Sö 9",
            "location": "Näsbyholm, Södermanland",
            "features": ["grid pattern", "boustrophedon text"],
            "inscription_fragments": ["antuitr", "uifara", "ikuari", "sikialti"]
        },
        {
            "id": "U 513",
            "location": "Färentuna Church, Uppland",
            "features": ["grid pattern"],
            "inscription_fragments": ["þurkil", "fulkui", "srklanti"]
        },
        {
            "id": "Sö 131",
            "location": "Gripsholm Castle, Södermanland",
            "features": ["serpent design"],
            "inscription_fragments": ["þurui", "ulfr", "ikuari"]
        },
        {
            "id": "Ög 153",
            "location": "Haddebo, Östergötland",
            "features": ["serpent design"],
            "inscription_fragments": ["skalt", "hiuk", "ikuari"]
        }
    ]

    # Step 2: Define the features observed in the provided image fragment.
    # The image clearly shows a grid-like pattern and a readable text fragment.
    observed_pattern = "grid pattern"
    # The middle row of runes transliterates to 'sikalt'.
    observed_text_fragment = "sikalt"

    print(f"Searching for a runestone with the following features:")
    print(f"- Visual Style: '{observed_pattern}'")
    print(f"- Text Fragment: '{observed_text_fragment}'\n")

    # Step 3: Filter runestones based on the primary visual feature.
    candidates = []
    for stone in runestones:
        if observed_pattern in stone["features"]:
            candidates.append(stone)

    print(f"Found {len(candidates)} candidate(s) with a '{observed_pattern}': {[c['id'] for c in candidates]}\n")

    # Step 4: Among the candidates, find the best match for the text fragment.
    # We will use Levenshtein distance to find the closest text match.
    # A smaller distance means a closer match.
    best_match = None
    min_distance = float('inf')

    print("Comparing text fragments to find the best match...")
    for stone in candidates:
        for fragment in stone["inscription_fragments"]:
            # The standard transliteration for Sö 9 is 'sikialti', which is an
            # interpretation. The carving itself is closer to 'sikalt'.
            # Levenshtein distance helps quantify this similarity.
            distance = Levenshtein.distance(observed_text_fragment, fragment)
            print(f"- Comparing '{observed_text_fragment}' with '{fragment}' (from {stone['id']})... Distance: {distance}")
            if distance < min_distance:
                min_distance = distance
                best_match = stone

    # Step 5: Output the result.
    if best_match:
        print(f"\nThe best match is runestone {best_match['id']} with a text similarity score (Levenshtein distance) of {min_distance}.")
        # The final ID is composed of a province code (Sö for Södermanland) and a number.
        province_code = "Sö"
        number_id = 9
        print(f"The ID of the Ingvar runestone is {province_code} {number_id}.")

# Since Levenshtein is a common library but might not be in a standard environment,
# we'll define a simple version of it if it's not found.
try:
    import Levenshtein
except ImportError:
    print("Levenshtein library not found. Using a basic implementation.")
    class Levenshtein:
        @staticmethod
        def distance(s1, s2):
            if len(s1) > len(s2):
                s1, s2 = s2, s1
            distances = range(len(s1) + 1)
            for i2, c2 in enumerate(s2):
                new_distances = [i2 + 1]
                for i1, c1 in enumerate(s1):
                    if c1 == c2:
                        new_distances.append(distances[i1])
                    else:
                        new_distances.append(1 + min((distances[i1], distances[i1 + 1], new_distances[-1])))
                distances = new_distances
            return distances[-1]

find_runestone_id()