def solve():
    """
    Solves the two parts of the user's question based on the reasoning above.
    """

    # Answer to Part A: Is every boopable function Turing-computable?
    answer_A = "Yes"

    # Answer to Part B: Does there exist a program L of length <= 2^32+50 such that,
    # for all x in [0, 2^32), L on (R,A,W) maps x to x and halts in <= x+50 many steps?
    answer_B = "No"
    
    print(f"{answer_A} {answer_B}")

solve()