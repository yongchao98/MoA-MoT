def solve():
    """
    Calculates the complexity for the two regimes.
    """
    # Regime 1: N = 2^sqrt(L)
    # The optimal strategy has complexity Theta(N log N).
    # sqrt(N^a (log N)^b (log log N)^c) = N log N
    # N^a (log N)^b (log log N)^c = N^2 (log N)^2
    # This gives (a,b,c) = (2,2,0).
    regime1_complexity = (2, 2, 0)

    # Regime 2: N = 2^((log L)^2)
    # The optimal strategy also has complexity Theta(N log N).
    # The reasoning is similar to regime 1.
    # This gives (a,b,c) = (2,2,0).
    regime2_complexity = (2, 2, 0)
    
    # Format the output as requested in the problem description example.
    # The format "(c1),(c2)" means we format the two tuples.
    r1_str = f"({regime1_complexity[0]},{regime1_complexity[1]},{regime1_complexity[2]})"
    r2_str = f"({regime2_complexity[0]},{regime2_complexity[1]},{regime2_complexity[2]})"
    
    print(f"{r1_str},{r2_str}")

solve()