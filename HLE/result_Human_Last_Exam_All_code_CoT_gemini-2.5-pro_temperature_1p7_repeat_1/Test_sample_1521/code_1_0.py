def solve_markov_chain_questions():
  """
  Solves the two theoretical questions about Markov chains.
  
  Question 1: Given an irreducible Markov chain and a function h with certain properties
  (harmonic outside a finite set A, zero on A, tending to infinity), can one conclude
  if the chain is recurrent or transient?
  
  Analysis 1: The existence of such a function h implies that the probability of never
  hitting the finite set A is positive. This is a defining characteristic of a transient
  chain, as a recurrent chain would visit any finite set with probability 1.
  So, the first answer is 't' for transient.

  Question 2: Considering a new chain which is a Doob's h-transform of the first,
  q(x,y) = p(x,y) * h(y)/h(x), is this new chain recurrent or transient?

  Analysis 2: The h-transform results in a new chain for which 1/h(x) is a non-negative
  harmonic function. Since h(x) is not constant (as it tends to infinity), 1/h(x) is also
  not constant. The existence of a non-constant non-negative harmonic function implies
  that the chain must be transient.
  So, the second answer is also 't'.

  The final answer is presented as a tuple (first_answer, second_answer).
  """
  
  # Answer to the first question ('t' for transient)
  first_answer = 't'
  
  # Answer to the second question ('t' for transient)
  second_answer = 't'
  
  # The problem asks for the output in the format (first answer, second answer)
  final_answer = (first_answer, second_answer)
  
  print(f"{final_answer}")

solve_markov_chain_questions()