def solve_pisa_curse(graduation_score):
  """
  Explains the superstitious remedy for a student who climbed the Leaning Tower of Pisa before graduating.
  
  Args:
    graduation_score: The student's final graduation score.
  """
  
  print("According to a Pisan superstition, a student who climbs the Leaning Tower before graduating will not graduate.")
  print("However, there is a way to reverse the bad luck!")
  print("\nThe remedy involves another building in the Piazza dei Miracoli: the Baptistery of St. John.")
  print("\nTo break the curse, after receiving your final grade, you must circle the Baptistery one time for each point of your graduation score.")
  print("\nSo, for a student with a final score of " + str(graduation_score) + ", the equation is:")
  print(f"1 trip around the Baptistery * {graduation_score} points = {graduation_score} circles to break the curse.")

# Let's assume a student graduated with a top score of 110.
hypothetical_score = 110
solve_pisa_curse(hypothetical_score)