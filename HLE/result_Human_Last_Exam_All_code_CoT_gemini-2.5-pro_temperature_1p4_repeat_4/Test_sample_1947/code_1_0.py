# Coefficients for the number of closed tree-like walks of length 6.
# c1 for e (number of edges, P_2)
# c2 for k (number of triangles, K_3)
# c3 for p (number of paths of length 3, P_4)
# c4 for sum(deg(v) choose 2) (number of paths of length 2, P_3)
# c5 for sum(deg(v) choose 3) (number of star graphs, K_{1,3})

c1 = 2
c2 = 0
c3 = 6
c4 = 12
c5 = 12

# Print the coefficients in the specified order.
print(f"{c1}")
print(f"{c2}")
print(f"{c3}")
print(f"{c4}")
print(f"{c5}")