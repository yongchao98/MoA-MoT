def solve():
    """
    This function prints the answers to the user's questions.
    """
    # Answer to Part A: Is every boopable function Turing-computable?
    answer_a = "Yes"

    # Answer to Part B: Does there exist a program L with the specified properties?
    answer_b = "Yes"
    
    # The reasoning for Part B involves satisfying the time constraint equation.
    # Total_Steps = Setup_Steps + Boop_Steps + Halt_Step
    # Setup_Steps + x + 1 <= x + 50
    # Setup_Steps <= 49
    # A self-modifying program can achieve this with a small number of setup steps,
    # for example, around 15. The final check is that 15 is indeed less than or equal to 49.

    print(f"{answer_a} {answer_b}")

solve()