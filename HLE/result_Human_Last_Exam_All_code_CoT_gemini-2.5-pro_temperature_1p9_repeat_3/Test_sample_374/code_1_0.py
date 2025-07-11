import math

def get_odd_part(n):
  """Computes the odd part of a positive integer n."""
  if n == 0:
    return 0
  # A simple way to remove all factors of 2
  while n % 2 == 0:
    n //= 2
  return n

def order_gl_str(n, q):
  """Returns a string representing the calculation of the order of GL(n, q)."""
  if n == 0:
    return "1"
  
  terms = []
  order = 1
  for i in range(n):
    term_val = q**n - q**i
    terms.append(f"({q**n} - {q**i})")
    order *= term_val

  return " * ".join(terms) + f" = {order}"

def solve():
  n_dim = 4
  q = 2

  print("The problem asks for the highest possible order of the inertial quotient E.")
  print("Theory tells us this is equivalent to finding the maximum possible order of an odd-order")
  print(f"subgroup of Aut(D) = GL({n_dim}, {q}).")
  print("-" * 60)

  print("Let H be a subgroup of GL(4, 2) of odd order. By Maschke's Theorem, its action")
  print("on the vector space V = (F_2)^4 is completely reducible. This means V decomposes")
  print("into a direct sum of irreducible submodules, and H must embed into a direct")
  print("product of smaller GL groups. We analyze cases based on these decompositions.")
  print("-" * 60)
  
  print("Case 1: The action of H is reducible.")
  print("This corresponds to partitions of the dimension 4 into at least two parts.")
  partitions_of_4 = [[3, 1], [2, 2], [2, 1, 1], [1, 1, 1, 1]]
  max_reducible_order = 0

  for p in partitions_of_4:
    print(f"\nAnalyzing partition {p}:")
    total_odd_order = 1
    for dim in p:
      gl_order = 1
      for i in range(dim):
        gl_order *= (q**dim - q**i)
      odd_part = get_odd_part(gl_order)
      total_odd_order *= odd_part
      print(f"  - Submodule of dimension {dim} -> GL({dim}, 2).")
      print(f"    Order of GL({dim}, 2) is {order_gl_str(dim, q)}.")
      print(f"    The odd part of this order is {odd_part}.")
    
    print(f"  The maximum odd order for partition {p} is the product of odd parts: {total_odd_order}.")
    if total_odd_order > max_reducible_order:
      max_reducible_order = total_odd_order

  print(f"\n=> The maximum possible order for H in the reducible case is {max_reducible_order}.")
  print("-" * 60)
  
  print("Case 2: The action of H is irreducible.")
  print("H is an irreducible subgroup of GL(4, 2).")
  gl4_order_val = 1
  for i in range(4): gl4_order_val *= (q**4 - q**i)
  gl4_odd_part = get_odd_part(gl4_order_val)
  
  print(f"The order of GL(4, 2) is {order_gl_str(4, q)}.")
  print(f"The odd part of this order is {gl4_odd_part}, which is 3^2 * 5 * 7.")
  print("The order of H must be a divisor of 315.")

  print("\nAnalyzing existence of irreducible subgroups with large odd order:")
  print("  - An irreducible abelian subgroup must be cyclic. The largest odd-order cyclic subgroup")
  print(f"    in GL(4, 2) has order dividing 2^4-1=15. So, the max abelian odd order is 15.")
  print("  - Arguments from group theory show that subgroups of order 35, 45, 63, etc., do not")
  print("    exist in GL(4, 2). For instance, a subgroup of order 45 must be a direct product")
  print("    C_5 x (C_3 x C_3), which cannot act irreducibly because it's non-cyclic and abelian.")
  
  max_irreducible_order = 15
  print(f"\n=> The maximum possible order for H in the irreducible case is {max_irreducible_order}.")
  print("-" * 60)
  
  final_answer = max(max_reducible_order, max_irreducible_order)
  print("Conclusion: The highest possible order for E is the maximum of the orders found in")
  print("the reducible and irreducible cases.")
  print(f"max({max_reducible_order}, {max_irreducible_order}) = {final_answer}.")

solve()