import math

def get_complexity_tuple(N_description, L_description):
    """
    This function symbolizes the logic to determine the complexity tuple (a,b,c).
    
    The analysis concludes that the optimal query complexity Q is O(N * min(log N, sqrt(L))).
    We need to find (a,b,c) such that Q^2 = O(N^a * (log N)^b * (log log N)^c).

    Let's analyze the term min(log N, sqrt(L)) for the two regimes.
    Note: For simplicity, log will denote log base 2.
    """
    
    # Regime 1: N = 2^sqrt(L)  =>  L = (log N)^2
    if N_description == "2^sqrt(L)":
        # In this regime, L = (log N)^2.
        # Thus, sqrt(L) = sqrt((log N)^2) = log N.
        # So, min(log N, sqrt(L)) = log N.
        # The complexity is Q = O(N log N).
        # Q^2 = O(N^2 * (log N)^2).
        # Comparing with N^a * (log N)^b * (log log N)^c, we get:
        a = 2
        b = 2
        c = 0
        return (a, b, c)

    # Regime 2: N = 2^((log L)^2)  => L = 2^sqrt(log N)
    if N_description == "2^((log L)^2)":
        # We need to compare log N with sqrt(L).
        # sqrt(L) = sqrt(2^sqrt(log N)) = 2^(0.5 * sqrt(log N)).
        # To compare log N and 2^(0.5 * sqrt(log N)), we can compare their logarithms.
        # log(log N) vs. log(2^(0.5 * sqrt(log N))) = 0.5 * sqrt(log N).
        # The function sqrt(x) grows faster than log(x). For large N,
        # 0.5 * sqrt(log N) > log(log N), which implies 2^(0.5 * sqrt(log N)) > log N.
        # Therefore, sqrt(L) > log N, and min(log N, sqrt(L)) = log N.
        # The complexity is Q = O(N log N).
        # Q^2 = O(N^2 * (log N)^2).
        # This gives the tuple:
        a = 2
        b = 2
        c = 0
        return (a, b, c)

# Calculate complexity for the first regime
regime1_tuple = get_complexity_tuple("2^sqrt(L)", "L")
# Calculate complexity for the second regime
regime2_tuple = get_complexity_tuple("2^((log L)^2)", "L")

# Format and print the final answer as requested
# The final equation is the tuple pair showing the complexity for each regime.
print(f"({regime1_tuple[0]},{regime1_tuple[1]},{regime1_tuple[2]}),({regime2_tuple[0]},{regime2_tuple[1]},{regime2_tuple[2]})")
