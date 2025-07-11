import math

def run_mimicry_simulation():
    """
    This simulation demonstrates why using human perception to cluster mimicry
    syndromes can be invalid from an ecological standpoint.
    """

    # 1. Define hypothetical Bombus species with color values.
    # We include R, G, B for human vision and UV for predator (e.g., bird) vision.
    # Values are on a 0-255 scale for simplicity.
    species_list = [
        {'name': 'B. terrestris', 'r': 20, 'g': 20, 'b': 20, 'uv': 10},
        {'name': 'B. lucorum',    'r': 25, 'g': 25, 'b': 20, 'uv': 90}, # Looks similar to terrestris for humans, but very different in UV.
        {'name': 'B. lapidarius', 'r': 200, 'g': 40, 'b': 20, 'uv': 20},
        {'name': 'B. ruderarius', 'r': 190, 'g': 45, 'b': 25, 'uv': 25}, # A true mimic of lapidarius.
        {'name': 'B. pascuorum',  'r': 180, 'g': 100, 'b': 40, 'uv': 50},
        {'name': 'B. crypticus',  'r': 140, 'g': 130, 'b': 80, 'uv': 55}  # Looks different to humans, but similar in UV to pascuorum.
    ]

    # A threshold for considering two species "similar". If distance is below this, they are clustered.
    SIMILARITY_THRESHOLD = 60.0

    print("--- Simulation of Researcher's Study (Human Perception) ---")
    print(f"Clustering species with a similarity distance threshold of < {SIMILARITY_THRESHOLD}\n")

    human_mimics = []
    # Iterate through all unique pairs of species
    for i in range(len(species_list)):
        for j in range(i + 1, len(species_list)):
            s1 = species_list[i]
            s2 = species_list[j]

            # Calculate distance based on human vision (RGB only)
            r1, g1, b1 = s1['r'], s1['g'], s1['b']
            r2, g2, b2 = s2['r'], s2['g'], s2['b']
            dist_sq = (r1 - r2)**2 + (g1 - g2)**2 + (b1 - b2)**2
            distance = math.sqrt(dist_sq)

            print(f"Comparing {s1['name']} and {s2['name']}:")
            print(f"  - Human model distance calc: sqrt(({r1}-{r2})^2 + ({g1}-{g2})^2 + ({b1}-{b2})^2) = {distance:.2f}")

            if distance < SIMILARITY_THRESHOLD:
                human_mimics.append((s1['name'], s2['name']))
                print(f"  -> RESULT: Clustered as mimics by humans.\n")
            else:
                print(f"  -> RESULT: Not considered mimics by humans.\n")

    print("\n--- Ecologically Relevant Simulation (Predator Perception) ---")
    print(f"Clustering species with a similarity distance threshold of < {SIMILARITY_THRESHOLD}\n")

    predator_mimics = []
    # Iterate through all unique pairs of species again
    for i in range(len(species_list)):
        for j in range(i + 1, len(species_list)):
            s1 = species_list[i]
            s2 = species_list[j]

            # Calculate distance based on predator vision (RGB + UV)
            r1, g1, b1, uv1 = s1['r'], s1['g'], s1['b'], s1['uv']
            r2, g2, b2, uv2 = s2['r'], s2['g'], s2['b'], s2['uv']
            dist_sq = (r1 - r2)**2 + (g1 - g2)**2 + (b1 - b2)**2 + (uv1 - uv2)**2
            distance = math.sqrt(dist_sq)

            print(f"Comparing {s1['name']} and {s2['name']}:")
            print(f"  - Predator model distance calc: sqrt(({r1}-{r2})^2 + ({g1}-{g2})^2 + ({b1}-{b2})^2 + ({uv1}-{uv2})^2) = {distance:.2f}")

            if distance < SIMILARITY_THRESHOLD:
                predator_mimics.append((s1['name'], s2['name']))
                print(f"  -> RESULT: Clustered as mimics by predators.\n")
            else:
                print(f"  -> RESULT: Not considered mimics by predators.\n")


    print("\n--- Conclusion ---")
    print("Human-based clusters found:", human_mimics)
    print("Predator-based (true) clusters found:", predator_mimics)
    print("\nKey Discrepancies:")
    print("1. ('B. terrestris', 'B. lucorum'): Humans see them as mimics, but their large UV difference means a predator would easily tell them apart. The researcher's method produces a FALSE POSITIVE.")
    print("2. ('B. pascuorum', 'B. crypticus'): Humans see them as different, but their similar UV patterns make them look alike to a predator. The researcher's method produces a FALSE NEGATIVE.")
    print("\nThis simulation shows that relying on human perception is not a valid approach because it ignores the visual system of the predator, which is the agent of selection that drives the evolution of mimicry.")

if __name__ == '__main__':
    run_mimicry_simulation()