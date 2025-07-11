import math

def calculate_s_curve_threshold():
    """
    Calculates the similarity threshold 's' for an LSH S-curve
    given r, b, and a target false positive rate.
    """
    r = 3
    b = 9
    fp_rate = 0.05

    print(f"Given parameters: r={r}, b={b}")
    print(f"The S-curve equation for the probability of collision is: P(s) = 1 - (1 - s^r)^b")
    print(f"We need to find the threshold 's' where the false positive rate is {fp_rate}.")
    print("\nSetting up the equation:")
    print(f"1 - (1 - s^{r})^{b} = {fp_rate}")
    
    # Step-by-step calculation
    val1 = 1 - fp_rate
    print(f"(1 - s^{r})^{b} = 1 - {fp_rate} = {val1}")

    val2 = math.pow(val1, 1/b)
    print(f"1 - s^{r} = {val1}^(1/{b})")

    val3 = 1 - val2
    print(f"s^{r} = 1 - {val1}^(1/{b})")

    s = math.pow(val3, 1/r)
    print(f"s = (1 - {val1}^(1/{b}))^(1/{r})")

    # Final result rounded
    s_rounded = round(s, 3)

    print(f"\nThe calculated threshold s is: {s_rounded}")

    return s_rounded

if __name__ == '__main__':
    result = calculate_s_curve_threshold()
    print(f"<<<{result}>>>")
