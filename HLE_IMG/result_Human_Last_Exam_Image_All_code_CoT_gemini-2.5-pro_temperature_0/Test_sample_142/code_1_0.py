import textwrap

def explain_beetle_observation():
    """
    Explains why observing the Rainbow Leaf Beetle in Germany is unlikely.
    """
    
    # Step 1: Identify the insect
    step1 = "Step 1: The insect in the image is the Rainbow Leaf Beetle (Chrysolina cerealis)."
    
    # Step 2: Determine its geographic range
    step2 = "Step 2: This beetle is native to Europe and parts of Asia, including Germany. This means options A (endemic to North America) and B (endemic to the tropics) are incorrect."
    
    # Step 3: Check its existence and conservation status
    step3 = "Step 3: The beetle is a real, living species, although it is rare in many areas. This rules out options D (it is not real) and E (it is extinct)."
    
    # Step 4: Evaluate the remaining options
    step4 = ("Step 4: We are left with options C and F. While the beetle is rare in Germany, the claim in F (not observed in over ten years) is likely false, as sporadic sightings are still recorded. Option C states its population has been reduced by over 76% in the last four decades.")
    
    # Step 5: Connect to real-world data
    step5_raw = ("Step 5: This figure aligns with a famous 2017 scientific study from Krefeld, Germany, which found a decline of more than 75% in the total biomass of flying insects in German nature reserves over just 27 years. This massive, documented decline in insects across Germany is the most significant reason why observing any specific insect, especially a rarer one like the Rainbow Leaf Beetle, has become very unlikely.")
    
    # Step 6: Conclude the answer
    step6 = "Step 6: Therefore, the widespread and severe decline in insect populations is the most accurate explanation provided."

    # Print the explanation
    print(step1)
    print(step2)
    print(step3)
    print(step4)
    print("\n".join(textwrap.wrap(step5_raw, width=80)))
    print(step6)
    
explain_beetle_observation()