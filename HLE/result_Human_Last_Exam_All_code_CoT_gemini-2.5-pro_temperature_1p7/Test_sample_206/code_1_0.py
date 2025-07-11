def solve_mtg_puzzle():
    """
    Solves the MTG puzzle by determining the optimal attack and resulting deaths.

    The logic is as follows:
    1.  Your life is 2, opponent's is 3. You must attack to win or you will likely lose next turn.
    2.  The goal is to deal 3 damage.
    3.  Your creatures: (1) Axegrinder Giant 6/4 Trample, (2) Centaur Courser 3/3, (3) Axebane Beast 2/5, (4) Wind Drake 2/2 Flying.
    4.  Opponent's creatures: (5) River Bear 4/3, (6) Grizzly Bears 2/2, (7) Rusted Sentinel 3/4, (8) Skywinder Drake 2/2 Flying (cannot block non-flyers).
    5.  Optimal Attack: The best course of action is an "all-out" attack to maximize pressure and damage, even if it's not lethal. This is the only move that puts the opponent in a position where a mistake could be lethal and it deals the maximum possible damage (2).
        - Attacking creatures: (1), (2), (3), (4).
    6.  Opponent's Optimal Block (to minimize damage taken):
        - (4) Wind Drake (2/2 F) is blocked by (8) Skywinder Drake (2/2 F). They trade and both die.
        - (1) Axegrinder Giant (6/4 T) is blocked by (7) Rusted Sentinel (3/4). This minimizes trample damage. The Giant deals 4 damage to kill the Sentinel, and 2 damage tramples over. The Sentinel dies.
        - (2) Centaur Courser (3/3) is blocked by (5) River Bear (4/3). They trade and both die.
        - (3) Axebane Beast (2/5) is blocked by (6) Grizzly Bears (2/2). The Beast survives, the Bears die.
    7.  Resulting deaths:
        - Your creatures: (2), (4).
        - Opponent's creatures: (5), (6), (7), (8).
    8.  The combined list of dying creatures is then sorted in increasing order.
    """
    
    attacking_creatures = [1, 2, 3, 4]
    
    # Creatures that die in the optimal block scenario
    dying_creatures = [2, 4, 5, 6, 7, 8]
    
    # Sort both lists to be sure, as per the format requirement
    attacking_creatures.sort()
    dying_creatures.sort()
    
    # Format the output string as per the example: (1), (2); (1), (6)
    attackers_str = ", ".join([f"({c})" for c in attacking_creatures])
    dying_str = ", ".join([f"({c})" for c in dying_creatures])
    
    final_answer = f"{attackers_str}; {dying_str}"
    
    print(final_answer)

solve_mtg_puzzle()
<<< (1), (2), (3), (4); (2), (4), (5), (6), (7), (8) >>>