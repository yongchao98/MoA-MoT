import math

def calculate_minimal_communication(p):
    """
    Calculates the minimal average one-way communication (C_1) required to simulate
    a singlet state, given an average of 'p' PR-boxes as a resource.

    Args:
        p (float): The probability of having a PR-box, between 0 and 1.

    Returns:
        float: The minimal average bits of one-way communication required.
    """
    if not (0 <= p <= 1):
        raise ValueError("The probability 'p' of having a PR-box must be between 0 and 1.")
    
    # The trade-off for minimal resources is given by the equation:
    # C_1 = (1 + sqrt(1 - p)) / 2
    min_C1 = (1 + math.sqrt(1 - p)) / 2
    return min_C1

def explain_and_demonstrate():
    """
    Explains the resource trade-off and demonstrates it with examples.
    """
    print("The minimal resources to classically simulate a singlet state involve a trade-off")
    print("between one-way communication (C_1) and PR-boxes (used with probability p).\n")
    print("The governing equation for the minimal resources is:")
    print("  C_1 = (1 + sqrt(1 - p)) / 2\n")

    # --- Case 1: No PR-boxes (p=0) ---
    p0 = 0.0
    c0 = calculate_minimal_communication(p0)
    print(f"Scenario 1: No PR-boxes are used (p = {p0}).")
    # Outputting the numbers in the equation for this case
    print(f"  Minimal C_1 = (1 + sqrt(1 - {p0})) / 2 = {c0:.2f} bits.")
    print("  This means exactly 1 bit of one-way communication is necessary and sufficient,")
    print("  a famous result by Toner and Bacon.\n")

    # --- Case 2: Intermediate case (e.g., p=0.75) ---
    p1 = 0.75
    c1 = calculate_minimal_communication(p1)
    print(f"Scenario 2: A PR-box is available 75% of the time (p = {p1}).")
    # Outputting the numbers in the equation for this case
    print(f"  Minimal C_1 = (1 + sqrt(1 - {p1})) / 2 = {c1:.4f} bits.")
    print(f"  The need for communication is reduced from 1 bit to {c1:.4f} bits on average.\n")

    # --- Case 3: Limit as p -> 1 ---
    p2_display = 1.0
    # Use a value very close to 1 for calculation to avoid domain errors if any
    p2_calc = 0.999999999 
    c2 = calculate_minimal_communication(p2_calc)
    print(f"Scenario 3: A PR-box is always available (limit as p -> {p2_display}).")
    # Outputting the numbers in the equation for this case
    print(f"  Minimal C_1 approaches (1 + sqrt(1 - {p2_display})) / 2 = {0.5:.2f} bits.")
    print("  This fundamental result shows that even with a perfect PR-box,")
    print("  0.5 bits of communication are still required for a perfect simulation.")

# Execute the explanation and demonstration
explain_and_demonstrate()