# The user wants a Yes/No answer to two questions.
# Question 1: Is there any a>0 real number for that floor(a^n) = n mod 2 for every n>0 integer?
# The analysis suggests that such a number 'a' can be constructed. The set of valid 'a' values is defined by an intersection of a nested sequence of non-empty intervals, which is non-empty.
answer1 = "Yes"

# Question 2: Is there any a>0 real number for that floor(a^n) = n mod 3 for every n>0 integer?
# The analysis via a constructive path led to a contradiction at n=5 for the most natural choice of floor values.
# This suggests that no such 'a' exists.
answer2 = "No"

# The final output should be the two answers, separated by a comma.
print(f"{answer1},{answer2}")
