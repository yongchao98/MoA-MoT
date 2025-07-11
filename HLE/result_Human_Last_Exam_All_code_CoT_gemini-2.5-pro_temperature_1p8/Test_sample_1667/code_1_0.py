def solve_duck_problem():
    """
    Solves the geometric probability problem about four ducks in a pond.

    This function explains the solution step-by-step using a symmetry argument
    and prints the final equation and result.
    """

    # Introduction to the problem's logic
    print("This problem can be solved elegantly using a symmetry argument rather than complex integration.")
    print("-" * 70)

    # Step 1 & 2: Define events and use symmetry
    print("Step 1: Define the Events and Use Symmetry\n")
    print("Let the four ducks be D1, D2, D3, and D4.")
    print("Let E_i be the event that duck Di is inside the circumcircle of the other three.")
    print("We want to find the probability of one of these events, let's say P(E4).\n")
    print("Because all ducks are placed under the same random conditions, the probability")
    print("of each event is the same. Let's call this probability 'p'.")
    print("So, p = P(E1) = P(E2) = P(E3) = P(E4).\n")
    print("-" * 70)

    # Step 3: State the key geometric property
    print("Step 2: A Key Geometric Property\n")
    print("For any four points in a plane (in a general position), it is a")
    print("mathematical fact that EXACTLY TWO of them lie inside the circumcircle")
    print("formed by the other three.\n")
    print("-" * 70)

    # Step 4: Apply Linearity of Expectation
    print("Step 3: Form an Equation Using Expected Values\n")
    print("Let N be the total number of 'in-circle' events that occur.")
    print("From the property above, we know that for any placement of ducks, N is always 2.")
    print("Therefore, the average or expected value of N, denoted E[N], must also be 2.\n")
    
    print("By the linearity of expectation, the expected value of N is also the sum of the")
    print("probabilities of the individual events:")
    print("E[N] = P(E1) + P(E2) + P(E3) + P(E4)")
    print("Using our symmetric probability 'p', this becomes:")
    print("E[N] = p + p + p + p = 4 * p\n")
    print("-" * 70)

    # Step 5: Solve the equation
    print("Step 4: Solve for the Probability 'p'\n")
    print("We now have two expressions for the expected value E[N]:")
    print("  1. E[N] = 2 (from the geometric property)")
    print("  2. E[N] = 4 * p (from the sum of probabilities)\n")
    
    print("Setting them equal gives us the final equation:")
    
    # Define numbers for the equation
    num_events = 4
    total_true_events = 2
    
    # Print the equation with each number explicitly
    print(f"  {num_events} * p = {total_true_events}\n")
    
    probability = total_true_events / num_events
    
    print("Solving for p, we get:")
    print(f"  p = {total_true_events} / {num_events} = {probability}\n")
    print("The probability that the fourth duck is within the circle formed by the first three is 0.5.")


# Run the function to display the solution
solve_duck_problem()

# The final answer
final_answer = 0.5