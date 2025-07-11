import math

# The problem is to find the limit of ln(1/p_n) / ln(n) as n -> infinity.
# Based on the derivation, the asymptotic behavior of p_n is C/n for some constant C.
# p_n ~ C / n
# So, 1/p_n ~ n / C
# ln(1/p_n) ~ ln(n / C) = ln(n) - ln(C)
# The expression becomes: (ln(n) - ln(C)) / ln(n) = 1 - ln(C)/ln(n)
# As n -> infinity, ln(n) -> infinity, so ln(C)/ln(n) -> 0.
# The limit is 1.

numerator_in_log = 1
denominator_in_log = 'p_n'
numerator_of_fraction = f"ln({numerator_in_log}/{denominator_in_log})"
denominator_of_fraction = "ln(n)"
final_answer = 1

print("The problem is to evaluate the limit:")
print(f"  lim ( {numerator_of_fraction} / {denominator_of_fraction} )")
print("  n->inf")
print("")
print("Based on the derivation, the asymptotic scaling of p_n is found to be:")
print("  p_n ~ C/n for some constant C")
print("")
print("Substituting this into the expression gives:")
print(f"  lim ( ln(n/C) / ln(n) ) = lim ( (ln(n) - ln(C)) / ln(n) ) = lim ( 1 - ln(C)/ln(n) )")
print("  n->inf")
print("")
print("The final result is:")
print(f"  {final_answer}")
