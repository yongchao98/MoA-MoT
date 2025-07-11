import sys

def solve_checkout_problem():
    """
    Calculates the expected checkout time for a friend based on the properties
    of memoryless distributions.
    """
    # Step 1: Define the known information from the problem.
    # We are told that the expected checkout time for a customer is 10 minutes.
    # This is the mean of the underlying probability distribution.
    expected_time_per_customer = 10

    # Step 2: Explain the underlying principle.
    # The problem describes a process where the chance of completion is constant over time.
    # This means the checkout duration follows an Exponential Distribution.
    # A key property of events governed by this distribution is that they are
    # independent. The duration of one customer's checkout does not affect the next.
    print("This problem describes a process whose duration follows an exponential distribution.")
    print("A key property is that each customer's checkout is an independent event.")
    print("\n--------------------------------------------------")

    # Step 3: Relate the principle to the specific question.
    print(f"The given expected checkout time for any customer is: {expected_time_per_customer} minutes.")
    print("Your friend is a new customer, so their checkout process is a new, independent event.")
    print("Therefore, their expected checkout time is drawn from the same distribution with the same mean.")
    print("--------------------------------------------------\n")

    # Step 4: State the final conclusion and the equation.
    expected_time_friend = expected_time_per_customer
    
    print("Final Answer Logic:")
    print("The expected time until your friend is done is equal to the average time for any customer.")
    
    # Final equation format as requested
    print("\nFinal Equation:")
    print(f"E[Friend's Time] = E[Any Customer's Time] = {expected_time_friend}")

solve_checkout_problem()
<<<C>>>