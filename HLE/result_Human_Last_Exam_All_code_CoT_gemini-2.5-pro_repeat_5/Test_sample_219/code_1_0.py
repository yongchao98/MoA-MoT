def analyze_system():
    """
    Analyzes the causal paths in a plant-pollinator system involving nectar caffeine.
    
    Variables:
    C: Nectar caffeine concentration
    F: Flower level foraging duration
    R: Pollinator retention
    Y: Total yield
    
    Paths:
    a: C -> F
    b: F -> Y
    c: C -> R
    d: R -> Y
    e: C -> Y (direct)
    """
    
    print("Analyzing the causal relationships in the plant-pollinator system...")
    print("="*60)
    
    # --- Analysis of each path ---
    
    # Path a: C -> F
    sign_a = '+'
    reasoning_a = (
        "Path 'a' (C -> F): Caffeine (C) is a stimulant. In nectar, it acts as a pharmacological reward for pollinators. "
        "Higher caffeine concentration encourages pollinators to spend more time foraging on a single flower (F). "
        "Therefore, the relationship is positive."
    )
    print(reasoning_a)
    print(f"Conclusion for 'a': {sign_a}\n")

    # Path b: F -> Y
    sign_b = '+'
    reasoning_b = (
        "Path 'b' (F -> Y): Longer flower-level foraging duration (F) typically leads to more effective pollination "
        "for that specific flower (i.e., more pollen is transferred). Better pollination increases the likelihood of "
        "successful fertilization and fruit set, thus increasing total yield (Y). Therefore, the relationship is positive."
    )
    print(reasoning_b)
    print(f"Conclusion for 'b': {sign_b}\n")

    # Path c: C -> R
    sign_c = '+'
    reasoning_c = (
        "Path 'c' (C -> R): Caffeine has been scientifically shown to enhance the memory of pollinators like bees. "
        "By remembering the location of the caffeinated nectar (C), pollinators are more likely to return to the same "
        "plants, increasing pollinator retention (R) and site fidelity. Therefore, the relationship is positive."
    )
    print(reasoning_c)
    print(f"Conclusion for 'c': {sign_c}\n")

    # Path d: R -> Y
    sign_d = '+'
    reasoning_d = (
        "Path 'd' (R -> Y): Higher pollinator retention (R) means a more consistent and dedicated workforce of pollinators "
        "for the plants. More frequent visits lead to higher rates of pollination and, consequently, a higher "
        "total yield (Y). Therefore, the relationship is positive."
    )
    print(reasoning_d)
    print(f"Conclusion for 'd': {sign_d}\n")

    # Path e: C -> Y
    sign_e = '+'
    reasoning_e = (
        "Path 'e' (C -> Y): This direct path represents effects not mediated by F or R. Plants evolve traits to maximize "
        "reproduction. While producing caffeine (C) has a metabolic cost, it can also deter nectar robbers or prevent "
        "microbial growth in the nectar, protecting the reward for legitimate pollinators. These effects would have a "
        "net positive impact on total yield (Y). Therefore, the most likely direct relationship is positive."
    )
    print(reasoning_e)
    print(f"Conclusion for 'e': {sign_e}\n")

    print("="*60)
    print("Final Result:")
    print("Based on the analysis, the most likely set of signs is:")
    print(f"a : {sign_a}, b: {sign_b}, c: {sign_c}, d: {sign_d}, e: {sign_e}")

# Run the analysis
analyze_system()