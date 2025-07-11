import math

def analyze_k_matching_complexity():
    """
    Analyzes the complexity of counting k-matchings to find the maximum k
    for which the problem is solvable in subcubic time.
    """
    print("This program determines the maximum integer k for which counting k-matchings")
    print("is solvable in subcubic time, O(|V|^(3-epsilon)), under standard assumptions.\n")

    print("The analysis relies on a conditional lower bound for the time complexity of")
    print("counting k-matchings, which is Omega(|V|^E), where the exponent E = k/2 + 1.\n")

    print("We are looking for the maximum integer k where the lower bound exponent E is less than 3.")
    
    target_exponent = 3

    # Check for k=3
    k_3 = 3
    exponent_3 = k_3 / 2 + 1
    print("-" * 40)
    print(f"Case 1: k = {k_3}")
    print(f"The theoretical lower bound for the time complexity exponent is:")
    print(f"    {k_3}/2 + 1 = {exponent_3}")
    print(f"Since the exponent {exponent_3} is less than {target_exponent}, a subcubic algorithm is considered possible.")
    print("Indeed, an O(|V|^omega) algorithm exists where omega is ~2.373.")

    # Check for k=4
    k_4 = 4
    exponent_4 = k_4 / 2 + 1
    print("-" * 40)
    print(f"Case 2: k = {k_4}")
    print(f"The theoretical lower bound for the time complexity exponent is:")
    print(f"    {k_4}/2 + 1 = {exponent_4}")
    print(f"Since the exponent {exponent_4} is not less than {target_exponent}, a subcubic algorithm")
    print("is not expected to exist under standard fine-grained complexity assumptions.")
    print("-" * 40)

    max_k = 3
    print(f"\nBased on this analysis, the maximum integer k for which counting k-matchings")
    print(f"can be done in subcubic time is {max_k}.")


if __name__ == '__main__':
    analyze_k_matching_complexity()
