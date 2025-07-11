def explain_difference():
    """
    Explains the main sporophyte difference between Grimmia montana and Grimmia curviseta.
    """
    title = "Distinguishing Feature between Grimmia montana and Grimmia curviseta"
    print(title)
    print("=" * len(title))

    distinction_point = "The primary difference lies in the seta (the capsule stalk) and the resulting orientation of the capsule."
    
    g_montana_desc = """
    Grimmia montana:
    - Seta: Short and STRAIGHT.
    - Capsule Orientation: Held upright (erect) and is typically immersed or only slightly emergent above the leaves.
    """
    
    g_curviseta_desc = """
    Grimmia curviseta:
    - Seta: Distinctly CURVED or arcuate (as the name 'curviseta' suggests).
    - Capsule Orientation: Nodding or hanging downwards, often nestled within the leaves due to the curved stalk.
    """
    
    summary = "\nIn summary, at a glance, you can distinguish them by looking for an UPRIGHT capsule (G. montana) versus a HANGING/NODDING capsule (G. curviseta)."

    print(distinction_point)
    print(g_montana_desc)
    print(g_curviseta_desc)
    print(summary)

if __name__ == "__main__":
    explain_difference()