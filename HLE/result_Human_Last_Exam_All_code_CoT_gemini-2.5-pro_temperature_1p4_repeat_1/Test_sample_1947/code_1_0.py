# The coefficients for the expression for the number of closed tree-like walks of length 6.
# c_1: coefficient for e (number of edges)
# c_2: coefficient for k (number of K_3 subgraphs)
# c_3: coefficient for p (number of P_4 subgraphs)
# c_4: coefficient for sum(deg(v) choose 2) (number of P_3 subgraphs)
# c_5: coefficient for sum(deg(v) choose 3) (number of K_1,3 subgraphs)

c1 = 2
c2 = 0
c3 = 6
c4 = 8
c5 = 6

# The final expression is:
# 2 * e + 0 * k + 6 * p + 8 * sum(deg(v) choose 2) + 6 * sum(deg(v) choose 3)
# We print the coefficients in the requested order.
print(f"c_1 = {c1}")
print(f"c_2 = {c2}")
print(f"c_3 = {c3}")
print(f"c_4 = {c4}")
print(f"c_5 = {c5}")
print(f"The final equation is:")
print(f"{c1} * e + {c2} * k + {c3} * p + {c4} * sum(deg(v) choose 2) + {c5} * sum(deg(v) choose 3)")
