import sys

# This script is designed to run in Python 3.
if sys.version_info[0] < 3:
    sys.exit("This script requires Python 3.")

def solve_freds_lie_puzzle():
    """
    Solves the logic puzzle by analyzing each statement.
    """
    print("Analyzing the logic puzzle step-by-step:\n")

    # Define statements and rules
    s1 = "My name is Fred; and if yesterday was after tomorrow, it would be Wednesday."
    s3 = "The total number of my friend in my class is even."
    s4 = "The number of the males is equal to number of the females in my class."
    s5 = "My best friend is older than me. He is 18 and I am 19."
    
    rules = [
        "1. Fred tells the truth on only one day of the week.",
        "2. Exactly one of the four statements is true.",
        "3. Fred said the one true statement on his truth day."
    ]

    print("--- The Rules ---")
    for rule in rules:
        print(rule)
    print("\n--- Statement Analysis ---")

    # Analyze S5
    print("Step 1: Analyzing S5 -> '{}'".format(s5))
    print("The statement claims Fred is younger than his friend, but gives the ages as Fred: 19, Friend: 18.")
    print("Since 19 > 18, the statement contradicts itself. Therefore, S5 is FALSE.")
    print("-" * 20)

    # Analyze S4 and S3
    print("Step 2: Analyzing S4 -> '{}' and S3 -> '{}'".format(s4, s3))
    print("If S4 were true, the number of males would equal the number of females.")
    print("In this case, the total number of friends (males + females) would be an even number, making S3 automatically true.")
    print("This would violate Rule 2 (only one statement is true). Therefore, S4 must be FALSE.")
    print("-" * 20)

    # Analyze S1
    print("Step 3: Analyzing S1 -> '{}'".format(s1))
    print("This is a conditional 'if-then' statement. The premise is 'if yesterday was after tomorrow'.")
    print("This premise is a logical impossibility, so it is always false.")
    print("In logic, any conditional statement with a false premise is considered a 'vacuous truth'. Therefore, S1 is TRUE.")
    print("-" * 20)
    
    # Synthesize results
    print("Step 4: Synthesizing the results.")
    print("We have determined the truth values:")
    print("  - S1: TRUE (by logic)")
    print("  - S3: Must be FALSE (as only one statement can be true)")
    print("  - S4: FALSE (to prevent S3 from being true)")
    print("  - S5: FALSE (by contradiction)")
    print("This fits Rule 2, as S1 is the single true statement.")
    print("-" * 20)
    
    # Final conclusion
    print("Step 5: Determining the day.")
    print("According to Rule 3, Fred said the true statement (S1) on his truth day.")
    print("The statement S1 concludes with '...it would be Wednesday.'")
    print("Since Fred is telling the truth, the factual part of his statement must be accurate.")
    print("He is using the logical puzzle to state that the day is Wednesday.")
    
    final_day = "Wednesday"
    print(f"\nTherefore, the day Fred tells the truth is {final_day}.")

solve_freds_lie_puzzle()
<<<C>>>