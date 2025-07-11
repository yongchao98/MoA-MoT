def main():
    """
    This script calculates the answer to the problem and prints it in the specified format.
    """
    # For part (a), the answer is a symbolic expression derived from the analysis.
    answer_a_expr = "t*(t-1)/2"

    # For part (b), we are given t=5. We calculate the result using the formula.
    t = 5
    
    # Show each number in the final equation as requested.
    t_minus_1 = t - 1
    numerator = t * t_minus_1
    result_b = numerator // 2

    # Print the full answer in the required format.
    print("The final answer is derived as follows:")
    print(f"(a) The maximum number of roots as an expression in t is: {answer_a_expr}")
    print(f"(b) For t = {t}, the calculation is: {t}*({t_minus_1}) / 2 = {numerator} / 2 = {result_b}")
    print("\nFormatted Answer:")
    print(f"(a) {answer_a_expr}; (b) {result_b}")

if __name__ == "__main__":
    main()