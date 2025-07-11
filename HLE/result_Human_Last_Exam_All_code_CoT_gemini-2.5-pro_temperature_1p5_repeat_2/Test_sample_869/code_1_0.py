import math

def main():
    """
    Calculates the total and favorable distributions to find the probability.
    """
    # Calculate S = 25! / (5!)^5
    # Numerator is 25!
    # Denominator is (5!)^5
    s_numerator = math.factorial(25)
    s_denominator = math.factorial(5) ** 5
    S = s_numerator // s_denominator

    # Calculate F = 5! * (10! / (2!)^5)
    # 5! for permutations of winning types
    # 10! / (2!)^5 for distributing the remaining items
    f_permutations = math.factorial(5)
    f_distributions_numerator = math.factorial(10)
    f_distributions_denominator = math.factorial(2) ** 5
    f_distributions = f_distributions_numerator // f_distributions_denominator
    F = f_permutations * f_distributions

    print("Step-by-step calculation:")
    print(f"Total number of ways to distribute the items, S = 25! / (5!)^5")
    print(f"S = {s_numerator} / {s_denominator}")
    print(f"S = {S}\n")

    print("Number of favorable distributions, F = 5! * (10! / (2!)^5)")
    print(f"F = {f_permutations} * ({f_distributions_numerator} / {f_distributions_denominator})")
    print(f"F = {f_permutations} * {f_distributions}")
    print(f"F = {F}\n")
    
    # The probability is P = F / S. We will output the equation with the calculated numbers.
    print("The probability P is F / S:")
    print(f"P = {F} / {S}")

if __name__ == "__main__":
    main()
