import math

def get_answer():
    """
    This function provides the solution to the user's question.
    """
    
    # (a) True or false: If G is not a union of K_2's, then ex(n; G, K_{1,t}-ind) = Theta(n).
    answer_a = "True"

    # (b) True or false: ex(n; sK_2, K_{1,t}-ind) is independent of n.
    answer_b = "True"
    
    # (c) For G ~ sK_2, express the upper bound for ex(n; sK_2, K_{1,t}-ind) in terms of s and t.
    # The derivation shows the upper bound is the maximum of two quantities.
    # The two quantities correspond to two types of extremal graphs.
    # The first quantity is the number of edges in K_{2s-1}, which is binom(2s-1, 2).
    # This can be written as (s-1)*(2*s-1).
    # The second quantity is derived from the G_{s,n} graph, which is valid for n < s+t-1.
    # Its maximum number of edges is at n=s+t-2, and is binom(s-1, 2) + (s-1)*(t-1).
    # This can be written as (s-1)*(s-2)/2 + (s-1)*(t-1).
    term1 = "(s-1)*(2*s-1)"
    term2 = "(s-1)*(s-2)/2 + (s-1)*(t-1)"
    answer_c = f"max({term1}, {term2})"

    final_answer = f"(a) {answer_a}; (b) {answer_b}; (c) {answer_c}"
    print(final_answer)

get_answer()