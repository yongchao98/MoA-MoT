import sys

def solve_mimicry_puzzle():
    """
    This function solves the insect mimicry puzzle by identifying and pairing
    mimic insects with the insects whose leaf damage they imitate.
    """

    # Panel A (Tortoise Beetle) mimics the linear scrape damage it creates itself.
    # The first letter is the mimic, the second is the source of the damage.
    pair1_mimic = 'A'
    pair1_causer = 'A'

    # Panel C (Moth) mimics the blotchy leaf damage caused by a mining larva (Panel B).
    pair2_mimic = 'C'
    pair2_causer = 'B'

    # Panel E (Leaf Insect) mimics the chewing damage caused by an insect like the katydid (Panel F).
    pair3_mimic = 'E'
    pair3_causer = 'F'

    # Format the pairs as requested.
    # The final answer consists of three pairs, showing each match.
    final_answer = f"{pair1_mimic}{pair1_causer}, {pair2_mimic}{pair2_causer}, {pair3_mimic}{pair3_causer}"
    
    print(final_answer)

solve_mimicry_puzzle()