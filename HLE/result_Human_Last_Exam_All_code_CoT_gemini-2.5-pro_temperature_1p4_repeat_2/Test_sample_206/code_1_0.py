def solve_mtg_puzzle():
    """
    This function determines the optimal attack and the outcome for the given MTG scenario.
    
    Analysis:
    The player is at 2 life, and the opponent is at 3 life. The opponent has enough creatures
    to block and prevent a lethal attack this turn. The best play is an "alpha strike"
    (attacking with all creatures) to force the opponent to lose their entire board,
    setting up a win on the next turn.

    Optimal Attack: Attack with all creatures (1, 2, 3, 4).

    Opponent's Optimal Blocks (to survive):
    - Attacker (4) Wind Drake (2/2 flying) is blocked by (8) Skywinder Drake (2/2 flying). Both die.
    - Attacker (2) Centaur Courser (3/3) is blocked by (5) River Bear (3/3). Both die.
    - Attacker (1) Axegrinder Giant (6/4) is blocked by (7) Rusted Sentinel (3/4). The Sentinel dies.
    - Attacker (3) Axebane Beast (2/5) is blocked by (6) Grizzly Bears (2/2). The Grizzly Bears die.
    """
    
    # List of creatures that should attack, in increasing order
    attacking_creatures = [1, 2, 3, 4]
    
    # List of all creatures that will die in combat, in increasing order
    # Player's creatures: (2) Centaur Courser, (4) Wind Drake
    # Opponent's creatures: (5) River Bear, (6) Grizzly Bears, (7) Rusted Sentinel, (8) Skywinder Drake
    dying_creatures = [2, 4, 5, 6, 7, 8]
    
    # Format the output string as per the user's request, e.g., (1), (2); (1), (6)
    attackers_str = ", ".join([f"({c})" for c in attacking_creatures])
    dying_str = ", ".join([f"({c})" for c in dying_creatures])
    
    final_answer = f"{attackers_str}; {dying_str}"
    
    print(final_answer)

solve_mtg_puzzle()