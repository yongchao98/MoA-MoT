# Number of ways to return to identity in k moves
# These are known values from the study of the Rubik's cube group
N4 = 324
N5 = 0
N6 = 10044

# Number of choices for a single move
moves = 12

# |A1|: Solved after 4 moves. The first 4 moves must form a sequence that solves the cube (N4 ways).
# The last 2 moves can be anything.
A1 = N4 * moves**2

# |A2|: Solved after 5 moves. N(5) is 0 because the graph is bipartite (no odd-length paths back to start).
A2 = N5 * moves

# |A3|: Solved after 6 moves.
A3 = N6

# |A1 intersect A3|: Solved after 4 moves AND after 6 moves.
# S4=I, S6=I. S6 = M6*M5*S4 = M6*M5. So M6*M5=I.
# There are N4 ways for the first 4 moves.
# There are 12 ways for M5, and M6 is then fixed (M6=M5_inverse).
A1_intersect_A3 = N4 * moves

# Using the Principle of Inclusion-Exclusion, simplified because other intersections are empty:
# Total = |A1| + |A2| + |A3| - |A1 intersect A3|
total_permutations = A1 + A2 + A3 - A1_intersect_A3

# Let's print the calculation step-by-step as requested
# The formula is also equivalent to 132 * N4 + 12 * N5 + N6
calc_A1 = f"{N4} * {moves**2}"
calc_A2 = f"{N5} * {moves}"
calc_A3 = f"{N6}"
calc_A1_intersect_A3 = f"{N4} * {moves}"
final_calc_str = f"{A1} + {A2} + {A3} - {A1_intersect_A3}"

print(f"The number of permutations returning to the original configuration is calculated by:")
print(f"|A1| + |A2| + |A3| - |A1 intersect A3|")
print(f"= (N(4) * 12^2) + (N(5) * 12) + N(6) - (N(4) * 12)")
print(f"= ({calc_A1}) + ({calc_A2}) + ({calc_A3}) - ({calc_A1_intersect_A3})")
print(f"= {final_calc_str}")
print(f"= {total_permutations}")

# Another way to express the final calculation:
# Total = N4 * (moves**2 - moves) + N5 * moves + N6
# Total = N4 * (144 - 12) + N5 * 12 + N6
# Total = N4 * 132 + N5 * 12 + N6
final_calc_simplified = f"{N4} * 132 + {N5} * 12 + {N6}"
print("\nOr, in simplified form:")
print(f"132 * N(4) + 12 * N(5) + N(6)")
print(f"= {final_calc_simplified}")
print(f"= {132 * N4 + 12 * N5 + N6}")