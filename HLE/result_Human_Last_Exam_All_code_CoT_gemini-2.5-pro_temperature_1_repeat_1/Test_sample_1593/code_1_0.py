# The analysis shows that for both specified regimes, the optimal algorithm
# achieves a query complexity of Theta(N log N).
#
# Query Complexity Q = Theta(N * log(N))
#
# To express this in the format Theta(sqrt(N^a * (log N)^b * (log log N)^c)),
# we square Q:
# Q^2 = Theta((N * log N)^2) = Theta(N^2 * (log N)^2 * (log log N)^0)
#
# From this, we can determine the integer coefficients:
# a = 2
# b = 2
# c = 0
#
# This result holds for both regimes:
# 1. N = 2^sqrt(L)
# 2. N = 2^((log2 L)^2)
#
# The final answer is formatted as (a1,b1,c1),(a2,b2,c2).

regime1_complexity = (2, 2, 0)
regime2_complexity = (2, 2, 0)

# Printing each number in the final formatted string as requested.
final_answer_string = f"({regime1_complexity[0]},{regime1_complexity[1]},{regime1_complexity[2]}),({regime2_complexity[0]},{regime2_complexity[1]},{regime2_complexity[2]})"
print(final_answer_string)