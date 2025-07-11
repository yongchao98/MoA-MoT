import sys

def solve_group_theory_problem():
    """
    Solves the given group theory problem.

    The problem asks two questions about a group G with specific properties.
    (a) Does a unique minimal supergroup hat(G) exist for solving p-nonsingular systems?
    (b) What is the maximum possible derived length of hat(G)?
    """

    # Part (a): Existence and uniqueness of hat(G)
    # The properties described for hat(G) are those of a p-completion or
    # p-localization of G. For solvable groups like G, the existence and
    # uniqueness (up to isomorphism) of such a minimal completion is a
    # standard result in combinatorial group theory.
    answer_a = "Yes"

    # Part (b): Maximum possible derived length of hat(G)
    # The derived length of the p-completion hat(G) is the same as the derived
    # length of G, which is dl(G).
    # The group G has a subnormal series G = G_1 ... G_{n+1} = {1} of length n
    # with abelian factors. This implies G is solvable.
    # From the series, we can deduce that dl(G) <= n.
    # To show this is the maximum possible length, we can provide an example.
    # The group of (n+1)x(n+1) upper unitriangular matrices over the integers
    # has a derived length of n and satisfies the conditions on its factors.
    # Its p-completion also has a derived length of n.
    # Therefore, the maximum possible derived length is n. Since n is given as
    # a parameter of the problem, the answer is an expression.
    answer_b = "n"

    # The prompt asks for the numbers in the final equation.
    # The equation is G_{n+1} = {1}. The numbers are n+1 and 1.
    # Let's formulate the final answer string as requested.
    # We will output n+1 for demonstration purposes as per the unusual instruction,
    # however the derived length is n. Let's make the primary answer n as per the derivation.
    final_answer_string = f"(a) {answer_a}; (b) {answer_b}"
    print(final_answer_string)


# Execute the solution
solve_group_theory_problem()

# The final answer in the required format
# Based on the problem description, the expression for (b) is 'n'.
# The user asked for the answer in a specific format at the end.
final_answer_for_prompt = "(a) Yes; (b) n"
# This final block is for the answer submission as requested by the user prompt.
sys.stdout.write(f"\n<<<{final_answer_for_prompt}>>>")
