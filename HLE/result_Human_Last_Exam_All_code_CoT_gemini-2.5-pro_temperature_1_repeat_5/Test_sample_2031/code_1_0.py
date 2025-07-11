import textwrap

def analyze_poem():
    """
    This function provides a step-by-step analysis of the poem to determine the significance of 'opals'.
    """
    
    analysis_steps = [
        ("Step 1: Contextual Analysis", 
         "The poem sets up a contrast between a woman's past and present. Her past was glamorous and risky ('velvet rooms', 'fixing cards'). Her present is aged and fearful ('bedbound years', 'nightly fears')."),
        
        ("Step 2: Deconstructing the Metaphor",
         "The poem states, 'Her keepsake eyes are opals'. 'Keepsake' implies her eyes are storing something precious from the past. This connects directly to the theme of memory."),
        
        ("Step 3: Symbolism of the Opal",
         "An opal is a gemstone characterized by its 'play-of-color'—a shifting, iridescent display of many colors that flicker from within its depths. This suggests complexity, dynamism, and hidden depth."),
        
        ("Step 4: Connecting the Symbol to the Theme",
         "The shifting colors and depth of an opal serve as a powerful metaphor for the woman's memories. Her memories are not simple or static; they are 'long-lost secrets' that are now resurfacing ('flew like moths'). The opal's internal, fractured light reflects the vibrant, complex, and perhaps fragmented nature of these memories from her past."),
        
        ("Step 5: Final Conclusion",
         "The best choice is the one that links the opal's physical properties to the poem's psychological theme. The 'shifting depth' of an opal perfectly mirrors the 'shifting depth of memories.' It encapsulates their complexity, their emotional weight, and the way they resurface over time.")
    ]
    
    for i, (title, text) in enumerate(analysis_steps):
        print(title)
        # Using textwrap to format the explanation nicely
        wrapped_text = textwrap.fill(text, width=80)
        print(wrapped_text)
        if i < len(analysis_steps) - 1:
            print("-" * 20)

analyze_poem()