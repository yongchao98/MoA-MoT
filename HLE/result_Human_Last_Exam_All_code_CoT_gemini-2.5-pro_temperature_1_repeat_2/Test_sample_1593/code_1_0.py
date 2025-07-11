# This script calculates and prints the complexity tuples for the two regimes.

# The query complexity is expressed as Theta(sqrt(N^a * (log N)^b * (log log N)^c)).
# We need to find the integers (a, b, c) for each regime.

# --- Regime 1: N = 2^sqrt(L) ---
# Our analysis shows the query complexity for this regime is Theta(N * log N).
# To match this to the given format:
# N * log N = sqrt(N^2 * (log N)^2) = sqrt(N^2 * (log N)^2 * (log log N)^0)
# Therefore, the parameters are:
a1 = 2
b1 = 2
c1 = 0
regime1_tuple = f"({a1},{b1},{c1})"

# --- Regime 2: N = 2^((log L)^2) ---
# Our analysis shows the query complexity for this regime is also Theta(N * log N).
# The reasoning is identical to the first regime.
# N * log N = sqrt(N^2 * (log N)^2) = sqrt(N^2 * (log N)^2 * (log log N)^0)
# Therefore, the parameters are:
a2 = 2
b2 = 2
c2 = 0
regime2_tuple = f"({a2},{b2},{c2})"

# The final answer combines the two tuples.
final_answer = f"{regime1_tuple},{regime2_tuple}"

print(final_answer)
<<<({a1},{b1},{c1}),({a2},{b2},{c2})>>>