import sys

def solve_fred_riddle():
    """
    This function solves the logic puzzle about Fred's lying day by analyzing the given statements.
    """

    print("Step 1: Analyze the truth value of each statement.")
    print("--------------------------------------------------")
    
    print("\nS1: 'My name is Fred; and if yesterday was after tomorrow, it would be Wednesday.'")
    print("  - The first part, 'My name is Fred', is true as the speaker is Fred.")
    print("  - The second part, 'if yesterday was after tomorrow...', is a conditional statement. The premise 'yesterday was after tomorrow' is a logical impossibility, making it always false.")
    print("  - In logic, an 'if-then' statement with a false premise is always true.")
    print("  - Therefore, S1 as a whole (True AND True) is always TRUE.")

    print("\nS5: 'My best friend is older than me. He is 18 and I am 19.'")
    print("  - This statement claims the friend is older, but then gives ages that prove the friend is younger (18 < 19).")
    print("  - Since it contradicts itself, statement S5 is always FALSE.")
    
    print("\nS4: 'The number of the males is equal to number of the females in my class.'")
    print("S3: 'The total number of my friend in my class is even.'")
    print("  - Let M be the number of males and F be the number of females.")
    print("  - If S4 is true, then M = F. The total number of friends is M + F = 2*M, which is always an even number.")
    print("  - This means that if S4 is TRUE, S3 must also be TRUE.")

    print("\nStep 2: Apply the condition 'Only one statement is true'.")
    print("--------------------------------------------------")
    print("  - We know for a fact that S1 is always true.")
    print("  - If S1 is true, and only one statement can be true, then S1 MUST be the single true statement.")
    print("  - This forces all other statements to be false. So, S3 must be false, S4 must be false, and S5 is already false.")
    print("  - This is a consistent conclusion (If S4 is false, S3 can also be false).")
    print("  - The definitive truth values are: S1=True, S3=False, S4=False, S5=False.")

    print("\nStep 3: Apply the condition 'Fred said it on the day when he tells the truth'.")
    print("--------------------------------------------------")
    print("  - 'It' refers to the one true statement, which we've identified as S1.")
    print("  - So, Fred uttered statement S1.")
    print("  - He did this on a day when he tells the truth. This is consistent, as S1 is a true statement.")

    print("\nStep 4: Deduce the specific day.")
    print("--------------------------------------------------")
    print("  - We need to find the specific day Fred was telling the truth. The only clue must be within the statement he made, S1.")
    print("  - S1: 'My name is Fred; and if yesterday was after tomorrow, it would be Wednesday.'")
    print("  - The statement is logically true on any day, but it's a riddle. In riddles, no information is irrelevant.")
    print("  - The only specific day mentioned is Wednesday. The most logical conclusion in a riddle context is that the statement is about the day it is uttered.")
    print("  - Therefore, the true statement S1 is telling us that the day in question is Wednesday.")
    print("  - Since Fred spoke on a truth-telling day, Wednesday must be a day he tells the truth.")

    print("\nConclusion:")
    print("The day Fred says the truth is Wednesday.")
    
# Execute the solver
solve_fred_riddle()