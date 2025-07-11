import textwrap

def identify_theorist():
    """
    Analyzes a surrealist collage to identify its primary theoretical influence
    from the early 20th century.
    """

    # Step 1: Analyze the Artistic Style and Technique
    style = "Surrealist Collage"
    style_description = "This artistic style, which became prominent in the 1920s, involves creating new realities from fragmented and juxtaposed elements. It aims to unlock the power of the unconscious mind."
    
    print("--- Step 1: Identifying the Artistic Style ---")
    print(f"Style: {style}")
    print("\n".join(textwrap.wrap(style_description, 70)))
    print("-" * 50)

    # Step 2: Analyze the Content and Themes
    themes = {
        "Hybrid Figure": "A central figure with a human body and an animal (antelope) head, explicitly referenced by the text 'developing hybrids'.",
        "Dream Logic": "The illogical combination of disparate elements like a modern postcard, mythological text, wild landscapes, and an animal-headed person.",
        "Mythological Reference": "The name 'Persephone' is present, invoking the Greek myth of descent into the underworld, duality, and trauma.",
        "Civilization vs. Instinct": "Text from a 'World's Finest Resort Hotel' (representing order, luxury, and civilization) is contrasted with the wild, animalistic, and natural elements in the image.",
        "The Unconscious": "The overall composition suggests an internal psychological landscape or a dream-state, filled with symbolic meaning."
    }

    print("--- Step 2: Analyzing Key Themes in the Collage ---")
    for theme, desc in themes.items():
        print(f"\n* {theme}:")
        print("  " + "\n  ".join(textwrap.wrap(desc, 68)))
    print("-" * 50)
    
    # Step 3: Connect Themes to a Theorist
    analysis = """
    The Surrealist movement was explicitly and profoundly influenced by the theories of one particular thinker. The themes of dream analysis, the power of the unconscious, the repression of instinctual desires by civilization, and the interpretation of symbols (like the hybrid figure as 'condensation' and Persephone's story as a trauma narrative) are all central tenets of his work. While Carl Jung's work on archetypes is also relevant (especially regarding Persephone), the entire artistic movement of Surrealism is historically and philosophically rooted in the psychoanalytic framework developed by this theorist.
    """
    
    theorist = "Sigmund Freud"

    print("--- Step 3: Connecting Themes to an Early 20th-Century Theorist ---")
    print("\n".join(textwrap.wrap(analysis, 70)))
    print("-" * 50)

    # Step 4: Final Conclusion
    print("--- Final Answer ---")
    print("The theorist whose ideas on the unconscious, dreams, and civilization are most clearly drawn upon in this surrealist image is:")
    print(theorist)


# Run the analysis
if __name__ == "__main__":
    identify_theorist()