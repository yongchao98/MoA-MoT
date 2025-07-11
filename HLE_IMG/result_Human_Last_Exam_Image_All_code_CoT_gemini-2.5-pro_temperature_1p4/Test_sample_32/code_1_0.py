def explain_difference():
    """
    This function explains the key distinguishing sporophyte feature
    between Grimmia montana and Grimmia curviseta.
    """
    species1 = {
        "name": "Grimmia montana",
        "seta_length": "very short",
        "seta_curvature": "straight",
        "capsule_position": "immersed or emergent (sitting among the leaves)"
    }

    species2 = {
        "name": "Grimmia curviseta",
        "seta_length": "long",
        "seta_curvature": "distinctly curved (arcuate)",
        "capsule_position": "exserted (held well above the leaves)"
    }

    print("The main difference in the sporophyte that helps distinguish Grimmia montana and Grimmia curviseta at a glance is the seta (the stalk supporting the capsule):\n")
    print(f"- {species1['name']}: Has a {species1['seta_length']}, {species1['seta_curvature']} seta, making the capsule {species1['capsule_position']}.")
    print(f"- {species2['name']}: Has a {species2['seta_length']}, {species2['seta_curvature']} seta, which holds the capsule in an {species2['capsule_position']} position.")
    print("\nTherefore, the presence of a long, curved seta in G. curviseta versus a short, straight seta in G. montana is the key visual difference.")

explain_difference()