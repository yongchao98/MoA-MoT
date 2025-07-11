def reveal_profession():
    """
    This function provides the answer to the historical question about Andre Vernault.
    """
    profession = "diamond cutter"
    explanation = (
        "Andre Vernault was a diamond cutter from Antwerp. He was not a spy, but rather "
        "part of a mission for the Belgian government-in-exile to smuggle industrial "
        "diamonds out of Belgium to the UK. These diamonds were vital for the Allied "
        "war effort (used in precision tools for manufacturing). He concealed his "
        "profession to protect his mission, which led British counterintelligence "
        "(MI5) to mistakenly believe he was a German agent."
    )
    
    print(f"Andre Vernault's true profession, which he concealed from interrogators, was: {profession}.")
    print("\n--- Additional Context ---")
    print(explanation)

reveal_profession()