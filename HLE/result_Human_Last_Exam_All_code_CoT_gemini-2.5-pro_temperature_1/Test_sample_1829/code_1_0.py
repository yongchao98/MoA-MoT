def solve_puzzle():
    """
    This function explains the step-by-step solution to the logic puzzle.
    """
    print("Step 1: Understanding Fred's behavior.")
    print("The puzzle states Fred has 'the day' (singular) when he tells the truth.")
    print("This means Fred has exactly ONE truth-telling day and SIX lie-telling days per week.\n")

    print("Step 2: Analyzing the statements.")
    print("S5) 'My best friend is older than me. He is 18 and I am 19.'")
    print("This statement is a self-contradiction. If his friend is 18 and he is 19, his friend is younger. Therefore, S5 is always FALSE.\n")

    print("S4) 'The number of the males is equal to number of the females in my class.'")
    print("S3) 'The total number of my friend in my class is even.'")
    print("If S4 is true (males = females), then the total number of friends is 2 * males, which is always even. This means S3 must also be true.")
    print("So, we have a logical rule: If S4 is true, then S3 must be true.\n")

    print("S1) 'My name is Fred; and if yesterday was after tomorrow, it would be Wednesday.'")
    print("The condition 'if yesterday was after tomorrow' is a logical impossibility. In puzzles like this, a statement of the form 'if [impossible], then [claim]' is a convoluted way of asserting the claim.")
    print("Therefore, we interpret S1 as being equivalent to the simpler statement: 'Today is Wednesday'.")
    print("The truth of this statement depends on the day it is said.\n")

    print("Step 3: Combining the rules and finding the solution.")
    print("We are told that exactly ONE of the four statements is true.")
    print("Since S5 is always false, the single true statement must be S1, S3, or S4.\n")

    print("Let's consider the possibilities:")
    print(" - Possibility A: S4 is the true statement.")
    print("   If S4 is true, our logic rule says S3 must also be true. This would mean there are two true statements, which violates the 'only one is true' rule. So, this is impossible.\n")

    print(" - Possibility B: S3 is the true statement.")
    print("   If S3 is true, then S1 and S4 must be false. Fred would have said the true statement (S3) on his truth day. But this doesn't help us identify which day that is. It could be any day except Wednesday (since S1, 'Today is Wednesday', must be false). The puzzle would have no unique solution.\n")

    print(" - Possibility C: S1 is the true statement.")
    print("   If S1 ('Today is Wednesday') is the one true statement, this fits all the rules.")
    print("   1. Is the 'only one is true' rule satisfied? Yes. S1 is true, and S3, S4, and S5 are false. This is consistent with our S4->S3 rule.")
    print("   2. When was the statement made? Fred says his true statement on his truth day. For the statement 'Today is Wednesday' to be true, it must be said on Wednesday.")
    print("   3. This means Fred's single truth-telling day must be Wednesday.\n")

    print("Conclusion: This is the only scenario that results in a unique, logical answer.")
    print("The day Fred tells the truth is Wednesday.")

solve_puzzle()