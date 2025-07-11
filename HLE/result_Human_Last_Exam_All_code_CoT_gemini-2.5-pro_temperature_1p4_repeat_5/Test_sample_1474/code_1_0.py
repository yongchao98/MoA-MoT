def find_common_imagery():
    """
    Analyzes potential shared cinematic imagery between directors Fritz Lang and William Friedkin
    based on a list of choices.
    """
    choices = {
        'A': 'Aboriginal masks',
        'B': 'Magic wands',
        'C': 'The first ever cyborgs on screen',
        'D': 'Bugs',
        'E': 'None of the above'
    }

    print("Investigating shared cinematic elements between Fritz Lang and William Friedkin...\n")

    # --- Choice A Analysis ---
    print(f"Step 1: Analyzing Choice A: '{choices['A']}'")
    print("   - In William Friedkin's 'The Exorcist', a key object is an ancient statue/amulet of the demon Pazuzu, found in Iraq. This is Mesopotamian, not Australian Aboriginal.")
    print("   - This specific type of mask is not a known recurring motif in Fritz Lang's films.")
    print("   - Conclusion: Choice A is incorrect.\n")

    # --- Choice B Analysis ---
    print(f"Step 2: Analyzing Choice B: '{choices['B']}'")
    print("   - The works of both directors are known for their dark, serious, and often gritty tones (Expressionism, crime, horror).")
    print("   - Fairy-tale elements like magic wands do not thematically align with their major films.")
    print("   - Conclusion: Choice B is incorrect.\n")

    # --- Choice C Analysis ---
    print(f"Step 3: Analyzing Choice C: '{choices['C']}'")
    print("   - Fritz Lang's 'Metropolis' (1927) famously features the Maschinenmensch (Machine-Person), widely considered one of cinema's first and most iconic robots or cyborgs.")
    print("   - However, William Friedkin's filmography does not contain cyborgs. The question requires the element to be in BOTH directors' works.")
    print("   - Conclusion: Choice C is incorrect.\n")

    # --- Choice D Analysis ---
    print(f"Step 4: Analyzing Choice D: '{choices['D']}'")
    print("   - Fritz Lang directed one of his early successful serials called 'The Spiders' ('Die Spinnen', 1919-1920). The villains are a secret, spider-themed criminal organization.")
    print("   - William Friedkin directed the psychological horror film 'Bug' (2006), where the plot is driven by a paranoia about an insect infestation. He also used flies as an atmospheric omen of demonic presence in 'The Exorcist'.")
    print("   - Spiders (arachnids) and bugs/insects represent a clear thematic and imagistic link.")
    print("   - Conclusion: Choice D is correct.\n")

    print("============================================")
    print("Final Determination: Based on the analysis, the element that appears in both directors' oeuvres is 'Bugs'.")
    print("============================================")

find_common_imagery()