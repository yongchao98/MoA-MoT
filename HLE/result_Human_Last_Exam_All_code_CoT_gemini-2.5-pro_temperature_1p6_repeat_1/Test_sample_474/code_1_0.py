import collections

def cluster_species(species_list, observers):
    """
    Clusters species based on a list of perceived traits.
    For example, a human observer might only perceive 'visible_color',
    while a bird predator perceives both 'visible_color' and 'uv_pattern'.
    """
    clusters = collections.defaultdict(list)
    for species in species_list:
        # Create a tuple of trait values that the observer can see
        perceived_traits = tuple(species[observer] for observer in observers)
        # Add the species name to the cluster defined by these traits
        clusters[perceived_traits].append(species['name'])
    return clusters

def print_clusters(cluster_type, clusters):
    """Helper function to print the clustering results nicely."""
    print(f"--- {cluster_type} ---")
    print(f"This clustering is based on the perception of a {cluster_type.split()[0]}.\n")
    if not clusters:
        print("No clusters formed.")
        return

    cluster_id = 1
    for traits, members in clusters.items():
        # The 'equation' of the cluster is the combination of its members
        # based on shared traits. We will print each member.
        print(f"Cluster #{cluster_id} (Traits: {traits}):")
        for member in members:
            print(f"  - {member}")
        cluster_id += 1
    print("-" * 25 + "\n")


# Define hypothetical Bombus species with visible and UV properties
species = [
    {'name': 'Bombus confusus', 'visible_color': 'Yellow-Banded', 'uv_pattern': 'UV-Dull'},
    {'name': 'Bombus hortorum', 'visible_color': 'Yellow-Banded', 'uv_pattern': 'UV-Bright'},
    {'name': 'Bombus ruderatus', 'visible_color': 'Yellow-Banded', 'uv_pattern': 'UV-Dull'},
    {'name': 'Bombus pratorum', 'visible_color': 'Red-Tailed', 'uv_pattern': 'UV-Bright'},
]

# --- Simulation of the Researcher's Method ---
# The researcher uses humans, who can only perceive 'visible_color'.
human_observers = ['visible_color']
human_based_clusters = cluster_species(species, human_observers)
print_clusters("Human Perception Clusters", human_based_clusters)


# --- Simulation of an Ecologically Valid Method ---
# An ecological approach considers the predator (a bird), which perceives
# both 'visible_color' and 'uv_pattern'.
bird_predator_observers = ['visible_color', 'uv_pattern']
predator_based_clusters = cluster_species(species, bird_predator_observers)
print_clusters("Bird Predator Perception Clusters", predator_based_clusters)

print("Conclusion: The human-based method creates a large, inaccurate cluster of all 'Yellow-Banded' bees.")
print("The ecologically valid predator-based method correctly separates them into two different mimicry rings based on their UV signals.")
