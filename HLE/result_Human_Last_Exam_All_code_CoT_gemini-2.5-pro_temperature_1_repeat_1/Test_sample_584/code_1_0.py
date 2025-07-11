import math

def calculate_lsh_threshold():
    """
    Calculates the LSH similarity threshold for a given S-curve and probability.
    """
    # Given parameters
    r = 3
    b = 9
    probability = 0.05

    # The S-curve equation is P(s) = 1 - (1 - s^r)^b
    # We need to find s such that P(s) = probability
    #
    # Derivation:
    # probability = 1 - (1 - s^r)^b
    # (1 - s^r)^b = 1 - probability
    # 1 - s^r = (1 - probability)^(1/b)
    # s^r = 1 - (1 - probability)^(1/b)
    # s = (1 - (1 - probability)^(1/b))^(1/r)

    # Calculate the threshold s
    s = (1 - (1 - probability)**(1/b))**(1/r)

    # Print the equation with the numbers plugged in
    print(f"The equation to solve for the similarity threshold 's' is:")
    print(f"s = (1 - (1 - {probability})^(1/{b}))^(1/{r})")

    # Print the final result rounded to three decimal places
    print(f"\nThe calculated similarity threshold is: {s:.3f}")

if __name__ == "__main__":
    calculate_lsh_threshold()