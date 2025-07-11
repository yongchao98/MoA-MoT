import itertools

def calculate_similarity(species1, species2, vision_model='human'):
    """
    Calculates a simple similarity score based on shared visual traits.
    """
    if vision_model == 'human':
        features1 = set(species1['human_visible_colors'])
        features2 = set(species2['human_visible_colors'])
    elif vision_model == 'predator':
        # Predator sees all colors, including UV.
        features1 = set(species1['human_visible_colors'] + [species1['uv_pattern']])
        features2 = set(species2['human_visible_colors'] + [species2['uv_pattern']])
    else:
        return 0.0

    intersection = len(features1.intersection(features2))
    union = len(features1.union(features2))
    
    # Jaccard Similarity
    return intersection / union if union > 0 else 0.0

def find_clusters(species_list, vision_model, similarity_threshold=0.6):
    """
    Forms clusters of species based on a similarity threshold.
    """
    # Create a graph where nodes are species and edges exist if similarity > threshold
    adj_list = {s['name']: [] for s in species_list}
    for s1, s2 in itertools.combinations(species_list, 2):
        similarity = calculate_similarity(s1, s2, vision_model)
        if similarity >= similarity_threshold:
            adj_list[s1['name']].append(s2['name'])
            adj_list[s2['name']].append(s1['name'])
            
    # Find connected components (our clusters)
    clusters = []
    visited = set()
    for species_name in adj_list:
        if species_name not in visited:
            cluster = []
            stack = [species_name]
            visited.add(species_name)
            while stack:
                current_node = stack.pop()
                cluster.append(current_node)
                for neighbor in adj_list[current_node]:
                    if neighbor not in visited:
                        visited.add(neighbor)
                        stack.append(neighbor)
            clusters.append(sorted(cluster))
            
    return sorted(clusters)

def main():
    """
    Main function to run the simulation and print results.
    """
    # Define a set of hypothetical Bombus species.
    # Note that Species B and D look identical to humans, but have different UV patterns.
    # Note that Species A and C look different to humans, but share a UV pattern.
    species = [
        {'name': 'Bombus A', 'human_visible_colors': ['Yellow', 'Black'], 'uv_pattern': 'UV-Bright'},
        {'name': 'Bombus B', 'human_visible_colors': ['Orange', 'Black'], 'uv_pattern': 'UV-Dull'},
        {'name': 'Bombus C', 'human_visible_colors': ['White', 'Black'],  'uv_pattern': 'UV-Bright'},
        {'name': 'Bombus D', 'human_visible_colors': ['Orange', 'Black'], 'uv_pattern': 'UV-Bright'},
    ]

    print("--- Simulating Bee Mimicry Clustering ---\n")
    print("Scenario: Four bee species with visible colors and invisible UV patterns.")
    for s in species:
        print(f"- {s['name']}: Visible={s['human_visible_colors']}, Predator-Only={s['uv_pattern']}")

    print("\n-------------------------------------------------")
    print("1. Clustering based on UNTRAINED HUMAN perception (Visible Colors Only)")
    print("-------------------------------------------------")
    human_clusters = find_clusters(species, vision_model='human')
    print("Human-perceived similarity would likely produce the following clusters:")
    for i, cluster in enumerate(human_clusters):
        print(f"  Cluster {i+1}: {cluster}")
    print("Note: Bombus B and D are clustered because they share the same visible colors.")
    print("Bombus A and C are separate because their visible colors are different.\n")

    print("-------------------------------------------------")
    print("2. Clustering based on PREDATOR perception (Visible + UV patterns)")
    print("-------------------------------------------------")
    predator_clusters = find_clusters(species, vision_model='predator')
    print("A predator that sees UV would perceive the 'true' ecological clusters:")
    for i, cluster in enumerate(predator_clusters):
        print(f"  Cluster {i+1}: {cluster}")
    print("Note: Bombus A, C, and D are clustered because they share the 'UV-Bright' signal.")
    print("Bombus B is now isolated because its 'UV-Dull' pattern makes it different.\n")
    
    print("=================================================")
    print("Conclusion: The clustering results are different.")
    print("Basing clusters on human perception leads to an incorrect understanding of the mimicry syndrome.")
    print("=================================================")

if __name__ == '__main__':
    main()
