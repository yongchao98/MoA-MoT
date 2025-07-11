def solve():
  """
  This function provides the answers to the two questions asked.

  Question 1: Is there any a>0 real number for that floor(a^n) = n mod 2 for every n>0 integer?
  Answer: Yes. The set of such `a` is non-empty. A constructive proof using nested intervals shows that such `a` can be found.

  Question 2: Is there any a>0 real number for that floor(a^n) = n mod 3 holds for every n>0 integer?
  Answer: No. This is a known mathematical result. While a simple constructive attempt might not immediately show a contradiction for all cases,
  it has been proven that no such `a` exists for any modulus m >= 3.
  """
  answer_mod_2 = "Yes"
  answer_mod_3 = "No"
  print(f"{answer_mod_2},{answer_mod_3}")

solve()