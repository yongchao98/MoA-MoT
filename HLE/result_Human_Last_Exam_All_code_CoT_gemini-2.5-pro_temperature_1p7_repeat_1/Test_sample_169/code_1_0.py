import random

def run_simulation_and_explanation():
    """
    This script simulates the effect of latitude on plant Volatile Organic Compound (VOC) diversity,
    based on the ecological principle of the "information arms-race" between plants and their enemies.
    """

    # --- Ecological Explanation ---
    print("Ecological Rationale:")
    print("1. Latitudinal Gradient: Biotic pressure from parasites and herbivores is highest in the tropics (low latitude) and decreases towards the poles (high latitude).")
    print("2. Alpha (α) Diversity: In high-pressure tropical environments, individual plants evolve more complex chemical arsenals (more VOC types) to defend against a wider array of enemies. This results in higher α-diversity.")
    print("3. Beta (β) Diversity: To prevent parasites from easily adapting and spreading, plant populations in the tropics evolve high chemical variation between individuals. This 'moving target' defense results in higher β-diversity.")
    print("Conclusion: As latitude decreases (e.g., from 60°N to 0°), both α-diversity and β-diversity are predicted to increase. This means the effect of latitude is negative for both.\n")
    
    # --- Simulation ---
    
    # Function to calculate Jaccard Dissimilarity (a measure of beta diversity)
    # Result: 0 means identical sets, 1 means completely different sets.
    def calculate_beta_diversity(set1, set2):
        # The number of shared items between the two sets
        intersection_size = len(set1.intersection(set2))
        # The total number of unique items in both sets combined
        union_size = len(set1.union(set2))
        
        if union_size == 0:
            return 0.0 # Dissimilarity is 0 if both plants have no VOCs
            
        jaccard_similarity = intersection_size / union_size
        jaccard_dissimilarity = 1 - jaccard_similarity
        return jaccard_dissimilarity, intersection_size, union_size

    # Define the sites at different latitudes
    latitudes = [60, 30, 0] # 60°N (Temperate), 30°N (Sub-tropical), 0° (Equatorial)
    
    # A hypothetical global pool of all possible VOCs a plant could make
    global_voc_pool = list(range(1000))

    print("--- Simulation Results ---")
    
    for lat in latitudes:
        # --- Alpha Diversity Calculation ---
        # We model α-diversity as the richness (number of unique VOCs) of a plant's chemical profile.
        # It increases as latitude decreases. At 60°N, richness is 20. At 0°N, richness is 50.
        alpha_richness = int(20 + (60 - lat) * 0.5)

        # --- Beta Diversity Calculation ---
        # We model β-diversity by simulating two plants and seeing how different they are.
        # The 'uniqueness_factor' controls how much chemical overlap there is between plants.
        # At 60°N, uniqueness is 0 (high overlap). At 0°N, uniqueness is 1 (low overlap).
        uniqueness_factor = (60 - lat) / 60.0

        # Simulate Plant 1's VOCs
        plant1_vocs = set(random.sample(global_voc_pool, alpha_richness))
        
        # Simulate Plant 2's VOCs based on its similarity to Plant 1
        num_shared_vocs = int(alpha_richness * (1 - uniqueness_factor))
        num_new_vocs = alpha_richness - num_shared_vocs
        
        # Take some VOCs from Plant 1's profile
        shared_vocs = set(random.sample(list(plant1_vocs), num_shared_vocs))
        
        # And add some new VOCs from the global pool that Plant 1 doesn't have
        available_new_vocs = [v for v in global_voc_pool if v not in plant1_vocs]
        new_vocs = set(random.sample(available_new_vocs, num_new_vocs))
        
        plant2_vocs = shared_vocs.union(new_vocs)

        # Calculate the β-diversity metric (Jaccard Dissimilarity)
        beta_div, intersection, union = calculate_beta_diversity(plant1_vocs, plant2_vocs)

        # --- Output Results for the site ---
        print(f"\nSite at {lat}° Latitude:")
        print(f"  - Simulated Alpha Diversity (VOC Richness per Plant): {alpha_richness}")
        print(f"  - Simulated Beta Diversity (Jaccard Dissimilarity between two plants):")
        # I am outputting the numbers used in the final beta diversity equation as requested.
        print(f"    - VOCs in Plant 1: {len(plant1_vocs)}")
        print(f"    - VOCs in Plant 2: {len(plant2_vocs)}")
        print(f"    - Shared VOCs (Intersection): {intersection}")
        print(f"    - Total Unique VOCs (Union): {union}")
        print(f"    - Equation: 1 - (Intersection / Union) = 1 - ({intersection} / {union})")
        print(f"    - Final Beta Diversity Score: {beta_div:.3f}")

if __name__ == "__main__":
    run_simulation_and_explanation()