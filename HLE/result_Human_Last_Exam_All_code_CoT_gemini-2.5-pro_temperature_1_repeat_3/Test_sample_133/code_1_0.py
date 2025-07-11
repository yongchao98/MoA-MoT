import sys

def solve():
    """
    This function analyzes the Connect 4 board to find the optimal move(s) for 'O'.

    The board state is as follows (rows 1-6 from top, columns a-g from left):
      a b c d e f g
    1 . . . . . . .
    2 . . . . . . .
    3 . . . . . . .
    4 . . . O O . .
    5 O . X O X X X
    6 X O O X X O X

    Player to move: O

    Analysis Steps:
    1.  A win in one move is not possible. There is no spot where 'O' can play to get four in a row immediately.

    2.  The fastest path to victory is to secure a win on the next turn. This is achieved by making a move that creates two or more simultaneous threats that the opponent ('X') cannot block at the same time.

    3.  Let's consider the key strategic positions for 'O':
        - A horizontal pair at d4 and e4.
        - A diagonal line of three at c6, d5, and e4.

    4.  Evaluate candidate move `f4`:
        - 'O' plays at `f4`.
        - Threat 1 (Horizontal): This creates a line of three 'O's: `d4, e4, f4`. This line now threatens to win at two open ends: `c4` and `g4`. Both `c4` and `g4` are valid moves for the next player, as the spots below them are occupied.
        - Threat 2 (Diagonal): The move at `f4` also makes the winning spot `f3` playable. The pre-existing diagonal `c6, d5, e4` is a line of three 'O's. Placing a piece at `f3` would win the game. Since `f4` is now occupied, `f3` is a valid next move.
        - Conclusion for `f4`: The move `f4` creates at least three winning threats for 'O's next turn (`c4`, `g4`, `f3`). 'X' can only block one. Therefore, `f4` guarantees a win on the next turn. It is an optimal move.

    5.  Evaluate candidate move `c4`:
        - 'O' plays at `c4`.
        - This creates a line of three 'O's: `c4, d4, e4`.
        - This line threatens to win at its open ends: `b4` and `f4`.
        - The threat at `f4` is valid, as `f5` is occupied.
        - However, the threat at `b4` is NOT valid as an immediate threat. The spot below it, `b5`, is empty. A piece cannot be played at `b4` until `b5` is filled.
        - Conclusion for `c4`: This move creates only one immediate threat (`f4`), which 'X' can block. It does not guarantee a win on the next turn and is therefore not optimal.

    6.  Final conclusion: The only move that guarantees the fastest possible win is `f4`.
    """
    optimal_moves = ["f4"]
    print(", ".join(optimal_moves))

solve()