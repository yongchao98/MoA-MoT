import textwrap

def explain_ballet_techniques():
    """
    Analyzes the pirouette preparation styles of major ballet schools
    to determine which pair uses allongé arms with bent knees from fourth position.
    """
    # Dictionary holding the analysis for each ballet school's technique.
    school_analysis = {
        "Paris Opera Ballet School (French School)":
            "The French style is renowned for its elegance, precision, and graceful, elongated lines. The term 'allongé' (meaning 'elongated') is French and is a key aesthetic. An allongé arm position in preparation for turns aligns perfectly with the school's emphasis on sophisticated, long lines.",
        "School of American Ballet (SAB - Balanchine Method)":
            "The Balanchine style is famous for a very specific pirouette preparation from a deep fourth position plié. A defining feature is the use of dynamic, un-rounded, 'allongé' arms to create speed and a sharp, expansive line, which is a hallmark of the technique.",
        "The Royal Ballet School (English School)":
            "Influenced heavily by the Cecchetti method, this school's pirouette preparations typically feature rounded arms. The focus is on creating a strong, contained frame to control the turn, rather than an extended, allongé line.",
        "Vaganova Academy (Russian School)":
            "The Vaganova method emphasizes strength and coordination, using a very structured preparation with powerful, rounded arms. This arm position is designed to generate maximum force and control from the torso and back.",
        "La Scala Theatre Ballet School (Italian - Cecchetti Method)":
            "Like the Royal Ballet, La Scala's technique is rooted in the Cecchetti method. It prioritizes strong, rounded arm positions for pirouette preparations to ensure stability and technical precision."
    }

    print("Analyzing Pirouette Preparation Styles by Ballet Institution:")
    print("-" * 70)
    for school, analysis in school_analysis.items():
        print(f"\nInstitution: {school}")
        # textwrap helps format the explanation nicely.
        wrapped_analysis = textwrap.fill(analysis, width=70)
        print(f"Technique Notes: {wrapped_analysis}")
    
    print("\n" + "-" * 70)
    print("Conclusion:")
    conclusion_text = (
        "Based on the analysis, the School of American Ballet (Balanchine) has a "
        "defining, well-known technique using allongé arms for pirouette preparation. "
        "The Paris Opera Ballet School (French School) also stylistically favors "
        "the elongated 'allongé' line in its movements. The other schools "
        "characteristically use more rounded arm positions for this specific preparation."
    )
    print(textwrap.fill(conclusion_text, width=70))
    print("\nTherefore, the correct pair is the Paris Opera Ballet School and the School of American Ballet.")

# Execute the function to print the explanation.
explain_ballet_techniques()