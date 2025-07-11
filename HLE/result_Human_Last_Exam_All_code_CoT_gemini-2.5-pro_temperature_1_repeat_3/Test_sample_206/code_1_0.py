def solve_mtg_puzzle():
    """
    This function determines the optimal attack and its outcome for the given MTG scenario.

    The plan is as follows:
    1. Analyze the board state, life totals, and creatures for both players.
    2. My life is at a critical 2, and the opponent has a 3/1 flyer ((8) Skywinder Drake), which is a lethal threat on their next turn. My only answer to it is my (4) Wind Drake.
    3. Determine if a lethal attack is possible. The opponent is at 3 life. However, they have enough blockers to prevent all damage from any possible combination of my attackers. A winning attack is not possible this turn.
    4. Since a winning attack is off the table, the optimal play is the one that best ensures my survival and improves my board position.
    5. The best strategic move is to eliminate the primary threat, the (8) Skywinder Drake. This can be achieved by attacking with my (4) Wind Drake.
    6. If I attack with (4), the opponent's optimal play is to block with their (8). Taking 2 damage and going to 1 life is too risky for them.
    7. In the resulting combat, my 2/2 Wind Drake and their 3/1 Skywinder Drake will trade, and both will die.
    8. Therefore, the creature to attack with is (4), and the creatures that will die are (4) and (8).
    """
    
    attacking_creatures = "(4)"
    dying_creatures = "(4), (8)"
    
    print(f"{attacking_creatures}; {dying_creatures}")

solve_mtg_puzzle()