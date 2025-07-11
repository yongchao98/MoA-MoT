def solve():
    """
    This function determines the classification for the growth rates f_1(n) and f_2(n).

    f_1(n) is the upper bound on the arboricity of the sampled graph G' when c=1.
    f_2(n) is the upper bound on the arboricity of the sampled graph G' when c=2.

    Based on the analysis:
    - For c=1, sampling with p_u = 1/d_u results in a graph G' where the arboricity is O(1).
      This corresponds to category 1.
    - For c=2, sampling with p_u = 1/d_u^2 results in a graph G' where the arboricity is also O(1),
      and is even more sparse than the c=1 case. This also corresponds to category 1.

    The final two-digit number is formed by concatenating the categories for f_1 and f_2.
    """
    f1_category = 1
    f2_category = 1
    
    # The problem asks to return a two-digit number.
    # We combine the two digits into a single number.
    final_answer = f1_category * 10 + f2_category
    
    # The final output needs to be in the format "<<<answer>>>"
    print(f"{final_answer}")

solve()