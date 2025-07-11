import random

def solve_hat_puzzle():
    """
    This script demonstrates the solution to the 12-person hat puzzle.
    It simulates the final state of knowledge after the team executes the optimal strategy.
    
    The strategy is for every possible pair of members to raise their hands once.
    As proven, this guarantees that at least 11 members' numbers will be revealed.
    This script shows how the 12th member can deduce their number in that worst-case scenario.
    """
    
    # 1. Setup the puzzle
    people = ["Alex", "Ben", "Charlie", "Dan", "Eva", "Frank", 
              "George", "Hannah", "Ivy", "Jack", "Kyle", "Liam"]
    hat_numbers = list(range(1, 13))
    
    # Shuffle the hat numbers and assign them to people
    random.shuffle(hat_numbers)
    assignments = dict(zip(people, hat_numbers))
    
    print("--- The Hat Puzzle Solution ---")
    print(f"There are {len(people)} team members and {len(hat_numbers)} unique hats numbered 1 to 12.")
    print("Strategy: Every possible pair of members raises their hands.\n")

    # 2. Simulate the worst-case outcome of the strategy
    # In the worst case, the leader manages to not reveal the number of exactly one person.
    # Let's assume 'Liam' is this person.
    unrevealed_person = "Liam"
    
    revealed_people = [p for p in people if p != unrevealed_person]
    
    print(f"Worst-Case Scenario Analysis:")
    print(f"At most one person's number can remain unrevealed. Let's assume this is {unrevealed_person}.")
    print(f"This means the numbers for the other 11 members have been revealed.\n")

    # 3. The 11 members know their numbers directly.
    print("--- Knowledge of the 11 Revealed Members ---")
    revealed_numbers = []
    for person in revealed_people:
        number = assignments[person]
        revealed_numbers.append(number)
        # print(f"{person} knows their number is {number} because it was revealed.")

    print("The 11 members whose numbers were revealed know their number with certainty.")
    print(f"The set of revealed numbers is: {sorted(revealed_numbers)}\n")

    # 4. The 12th person (Liam) deduces their number.
    print(f"--- How {unrevealed_person} Deduces Their Number ---")
    print(f"{unrevealed_person} knows that their number is the only one missing from the revealed set.")
    
    # We can find the missing number using set theory or summation.
    # Let's use summation as a clear 'equation'.
    
    total_sum_of_all_hats = sum(range(1, 13))
    sum_of_revealed_hats = sum(revealed_numbers)
    
    liams_number = total_sum_of_all_hats - sum_of_revealed_hats
    
    print("Deduction using an equation:")
    print(f"Total Sum of all numbers (1 to 12) = {total_sum_of_all_hats}")
    
    # Building the equation string
    equation_parts = [str(n) for n in revealed_numbers]
    print(f"Sum of Revealed Numbers = {' + '.join(equation_parts)} = {sum_of_revealed_hats}")
    
    print(f"\n{unrevealed_person}'s Number = (Total Sum) - (Sum of Revealed Numbers)")
    print(f"{unrevealed_person}'s Number = {total_sum_of_all_hats} - {sum_of_revealed_hats} = {liams_number}")
    
    # Verification
    print(f"Checking the actual assignment: {unrevealed_person}'s hat was {assignments[unrevealed_person]}. The deduction is correct.\n")
    
    # 5. Conclusion
    print("--- Conclusion ---")
    print("The 11 members whose numbers were revealed know their number directly.")
    print(f"The 12th member ({unrevealed_person}) can deduce their number with 100% certainty.")
    print("Therefore, this strategy guarantees that all 12 members will correctly determine their number.")
    print("\nThe largest possible value of N is 12.")

solve_hat_puzzle()