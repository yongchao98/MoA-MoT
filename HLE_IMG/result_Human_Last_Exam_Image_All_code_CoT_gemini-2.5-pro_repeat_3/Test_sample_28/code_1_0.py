import textwrap

def solve_entomology_question():
    """
    Identifies the insect and determines the most likely collection locality
    based on its known distribution and the context provided.
    """

    # Step 1: Identify the insect
    insect_name = "Lycorma delicatula"
    common_name = "Spotted Lanternfly"
    identification_notes = (
        "The insect in the image is the Spotted Lanternfly. Key features visible are the "
        "bright red hindwings with black spots, which are displayed when the insect is startled. "
        "This coloration is a hallmark of the species."
    )

    # Step 2: Determine Geographic Distribution
    native_range = "China, India, and Vietnam."
    invasive_range = (
        "The Spotted Lanternfly is a highly invasive species. It was first detected in the "
        "United States in Berks County, Pennsylvania, in 2014. Since then, it has spread "
        "extensively throughout the northeastern U.S. It is also invasive in South Korea and Japan."
    )

    # Step 3: Analyze the Answer Choices
    options = {
        'A': 'Philadelphia, Pennsylvania, USA',
        'B': 'Buffalo, New York, USA',
        'C': 'Miami, Florida, USA',
        'D': 'Thimphu, Bhutan',
        'E': 'Munich, Bavaria, Germany',
        'F': 'Luodong, Taiwan',
        'G': 'Las Vegas, Nevada, USA',
        'H': 'Jinan, Shandong Province, China',
        'I': 'Baltimore, Maryland, USA',
        'J': 'Hengchun, Taiwan'
    }

    analysis = {
        'A': "Very likely. Philadelphia is near the original US detection site and has had a major, well-documented infestation for years.",
        'B': "Plausible. The insect is established in Buffalo, but Philadelphia is more famously associated with the initial outbreak.",
        'C': "Unlikely. Not known to be established in Miami.",
        'D': "Unlikely. Not a primary location in its native or invasive range.",
        'E': "Very unlikely. Not established in Germany.",
        'F': "Plausible, but less likely than the US invasive range for a US entomologist.",
        'G': "Very unlikely. Unsuitable arid climate.",
        'H': "Plausible. This is within the insect's native range.",
        'I': "Likely. Baltimore is also heavily infested, similar to Philadelphia.",
        'J': "Plausible, same as F."
    }

    # Step 4: Synthesize and Conclude
    conclusion = (
        "Comparing the likely options (A, B, H, I), Philadelphia (A) stands out. "
        "It was the first major metropolitan area in the U.S. to be severely impacted "
        "and is often considered 'ground zero' for the American invasion. For an entomologist "
        "from the Smithsonian (in nearby Washington D.C.), a collecting trip to the "
        "epicenter of this significant biological invasion in June (when adults are active) "
        "is extremely probable. Therefore, Philadelphia is the most likely collection locality."
    )

    # Print the detailed reasoning
    print("--- Insect Identification and Analysis ---")
    print(f"Insect: {insect_name} ({common_name})")
    print("\nIdentification Notes:")
    print(textwrap.fill(identification_notes, 80))
    print("\nGeographic Distribution:")
    print(f"Native Range: {native_range}")
    print(f"Invasive Range (Summary): {textwrap.fill(invasive_range, 80)}")
    print("\n--- Analysis of Localities ---")
    for key, loc in options.items():
        print(f"{key}. {loc}: {analysis[key]}")
    print("\n--- Conclusion ---")
    print(textwrap.fill(conclusion, 80))
    print("\nFinal Answer Choice: A")

solve_entomology_question()