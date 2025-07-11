# The number of closed tree-like walks of length 6 in a simple graph X
# can be written as an expression of the form:
# c1*e + c2*k + c3*p + c4*Sum(deg(v) choose 2) + c5*Sum(deg(v) choose 3)
# Based on combinatorial counting, we determined the integer coefficients.

# Coefficient for e (number of edges)
c1 = 2

# Coefficient for k (number of K_3 subgraphs)
c2 = -36

# Coefficient for p (number of P_4 subgraphs)
c3 = 6

# Coefficient for Sum(deg(v) choose 2) (number of P_3 subgraphs + 3*k)
c4 = 12

# Coefficient for Sum(deg(v) choose 3) (number of K_{1,3} subgraphs)
c5 = 12

# The problem asks to write the coefficients in order.
# Here is the final list of coefficients.
coefficients = [c1, c2, c3, c4, c5]

print("The coefficients are:")
print(f"c_1 = {c1}")
print(f"c_2 = {c2}")
print(f"c_3 = {c3}")
print(f"c_4 = {c4}")
print(f"c_5 = {c5}")

# Final Answer in the requested format
final_answer_string = ", ".join(map(str, coefficients))
# The final answer in the requested format.
# print(f"<<<{final_answer_string}>>>")