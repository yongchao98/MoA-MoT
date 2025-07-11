import math

def calculate_lsh_threshold():
    """
    Calculates the LSH similarity threshold for a given false positive rate.
    
    The S-curve function is P(s) = 1 - (1 - s^r)^b, where P(s) is the
    probability of a pair with similarity s becoming a candidate pair.
    We need to find the threshold s for which this probability is 0.05.
    """
    r = 3
    b = 9
    false_positive_rate = 0.05

    # We need to solve the equation:
    # false_positive_rate = 1 - (1 - s^r)^b
    # Rearranging for s gives:
    # s = (1 - (1 - false_positive_rate)^(1/b))^(1/r)
    
    threshold = (1 - math.pow(1 - false_positive_rate, 1 / b)) ** (1 / r)

    # Print the equation with the given values
    print("The S-curve equation is: P(s) = 1 - (1 - s^r)^b")
    print("\nGiven the parameters:")
    print(f"r = {r}")
    print(f"b = {b}")
    print(f"False positive rate (P(s)) = {false_positive_rate}")

    print("\nWe solve for 's' in the following equation:")
    print(f"{false_positive_rate} = 1 - (1 - s^{r})^{b}")
    print("\nRearranging the equation to solve for s:")
    print(f"s = (1 - (1 - {false_positive_rate})^(1/{b}))^(1/{r})")

    # Print the final result
    print(f"\nThe calculated similarity threshold 's' is: {threshold:.3f}")

if __name__ == "__main__":
    calculate_lsh_threshold()