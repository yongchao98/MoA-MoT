# Based on the analysis, the query complexity in both regimes is determined to be Theta(N * log(N)).
# Here is a summary of the derivation:
#
# 1. Define two primary sorting strategies:
#    a) Direct Sort: Use H-queries on full strings to identify unique ones, then C-queries to sort them.
#       This has a complexity of Theta(N * log(N)) as the number of unique strings is close to N in both regimes.
#    b) Radix Sort: Break strings into chunks of size l and sort in L/l passes. The optimized cost for this
#       strategy is found to be Theta(N * L / log(N)).
#
# 2. Compare strategies for each regime:
#    - Regime 1 (N = 2^sqrt(L) => L = (log N)^2):
#      The Radix Sort complexity is Theta(N * (log N)^2 / log(N)) = Theta(N * log(N)).
#      This is the same as Direct Sort. Thus, the optimal complexity is Theta(N * log(N)).
#
#    - Regime 2 (N = 2^((log L)^2) => L = 2^sqrt(log N)):
#      The Radix Sort complexity is Theta(N * 2^sqrt(log N) / log(N)).
#      This is asymptotically larger than the Direct Sort complexity of Theta(N * log(N)).
#      Therefore, the optimal strategy is Direct Sort, and the complexity is Theta(N * log(N)).
#
# 3. Convert Theta(N * log(N)) to the (a,b,c) format.
#    We are given the complexity class Q = Theta(sqrt(N^a * (log N)^b * (log log N)^c)).
#    Our derived complexity is Q = N * log(N).
#    Squaring our complexity gives Q^2 = (N * log(N))^2 = N^2 * (log N)^2.
#    By comparing terms with the squared given class, we find the final equation for the exponents:
#    N^a * (log N)^b * (log log N)^c = N^2 * (log N)^2 * (log log N)^0
#
#    This implies:
#    a = 2
#    b = 2
#    c = 0

# Set the calculated values for (a,b,c) for both regimes.
regime1_a = 2
regime1_b = 2
regime1_c = 0

regime2_a = 2
regime2_b = 2
regime2_c = 0

# Format the final answer string as requested.
answer = f"({regime1_a},{regime1_b},{regime1_c}),({regime2_a},{regime2_b},{regime2_c})"

print(answer)

# The final answer in the required format.
print("<<<({a1},{b1},{c1}),({a2},{b2},{c2})>>>".format(a1=regime1_a, b1=regime1_b, c1=regime1_c, a2=regime2_a, b2=regime2_b, c2=regime2_c))