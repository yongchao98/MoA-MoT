import textwrap

def solve_chesterton_puzzle():
    """
    This function explains and solves the Father Brown literary puzzle.
    """
    
    explanation = textwrap.dedent("""
    The puzzle refers to a passage in the introductory chapter of G.K. Chesterton's "The Secret of Father Brown".

    In this chapter, Father Brown explains the source of his insight into human nature. He recounts a confession from a man whose marriage had deteriorated, leading his wife to commit a crime (forgery).

    Father Brown uses two place names as metaphors to describe the transformation of the couple's love:
    - The beginning of their marriage was idyllic and romantic, like a runaway wedding.
    - The end was filled with quarrels and bitterness.
    
    The two place names he uses both begin with paired consonants.
    """).strip()

    place1_name = "Gretna Green"
    place1_explanation = "This town in Scotland was historically a destination for young English couples seeking to marry without parental consent. It symbolizes the romantic, passionate beginning of a relationship. The paired consonant is 'Gr'."
    
    place2_name = "Stony Stratford"
    place2_explanation = "This is a town in England. A 'Stony Stratford row' is an old phrase for a bitter quarrel or a sad ending, contrasting with the happy beginning. The paired consonant is 'St'."
    
    print(explanation)
    print("\nThe two place names are:")
    print("-" * 25)
    print(f"1. {place1_name}")
    print(f"   - Explanation: {place1_explanation}")
    print(f"2. {place2_name}")
    print(f"   - Explanation: {place2_explanation}")

# Execute the function to print the solution
solve_chesterton_puzzle()
