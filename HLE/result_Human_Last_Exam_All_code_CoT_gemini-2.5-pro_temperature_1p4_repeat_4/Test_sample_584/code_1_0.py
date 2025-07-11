import math

def calculate_lsh_threshold():
    """
    Calculates the LSH similarity threshold for a given false positive rate.
    """
    # Given parameters
    r = 3
    b = 9
    false_positive_rate = 0.05

    # The S-curve formula is P(s) = 1 - (1 - s^r)^b
    # We solve for 's' where P(s) = false_positive_rate
    #
    # false_positive_rate = 1 - (1 - s^r)^b
    # 1 - false_positive_rate = (1 - s^r)^b
    # (1 - false_positive_rate)^(1/b) = 1 - s^r
    # s^r = 1 - (1 - false_positive_rate)^(1/b)
    # s = (1 - (1 - false_positive_rate)^(1/b))^(1/r)

    # Calculate the threshold 's'
    s_r = 1 - (1 - false_positive_rate)**(1/b)
    s = s_r**(1/r)

    # Print the explanation and the final equation with values
    print("The goal is to find the similarity threshold 's' where the S-curve's value equals the desired false positive rate.")
    print(f"Given parameters r = {r}, b = {b}, and a false positive rate = {false_positive_rate}.")
    print("\nWe are solving the following equation for 's':")
    print(f"{false_positive_rate} = 1 - (1 - s^{r})^{b}")

    # Print the final calculated threshold
    print(f"\nThe similarity threshold 's' should be approximately: {s:.3f}")

if __name__ == "__main__":
    calculate_lsh_threshold()