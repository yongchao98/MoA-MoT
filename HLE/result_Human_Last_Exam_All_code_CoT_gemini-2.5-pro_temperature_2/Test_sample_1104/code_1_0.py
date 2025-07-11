import math

def solve_election_proportionality():
    """
    Calculates the smallest preference profile sizes (s1, s2) based on the problem's
    proportionality constraints.
    """
    
    # --- Problem Definition ---
    # Committee size
    k = 100
    
    # The initial 6 ballots define a critical group of voters.
    # A(1) = {a,b,c,x}, A(2) = {a,b,c,y}, A(3) = {a,b,c,y},
    # A(4) = {a,b,c,z}, A(5) = {a,b,c,z}, A(6) = {a,b,c,z}.
    # This group unanimously approves of {a,b,c}.
    
    # The size of this critical group.
    l = 6
    
    # The condition is that voter 1 is unsatisfied, which means {a,b,c,x} has an
    # empty intersection with the winning committee W. This implies that the
    # set {a,b,c} also has an empty intersection with W.

    print("Step 1: Analyzing the conditions for Proportional Justified Representation (PJR).")
    print(f"A critical group of voters of size l = {l} unanimously agrees on the set {{a,b,c}}.")
    print("For voter 1 to be unsatisfied, none of these candidates can be in the committee.")
    print("PJR states that if a group is 'large enough', it is owed a representative. A group is large enough if l >= s1/k.")
    print("To satisfy PJR, our group must NOT be large enough, because it doesn't get a representative.")
    print("This leads to the following inequality:")
    print(f"l < s1 / k")
    print(f"Substituting the values: {l} < s1 / {k}")
    s1 = l * k + 1
    print(f"From this, we get s1 > {l*k}. The smallest integer value for s1 is {s1}.")
    print("-" * 20)

    print("Step 2: Analyzing the conditions for Extended Justified Representation (EJR).")
    print(f"We consider the same critical group of size l = {l}.")
    print("This group gets 0 representatives from their common set {{a,b,c}}.")
    print("EJR states that if l >= j * (s2/k) for some j>=1, the group is owed at least j representatives.")
    print("To satisfy EJR, our group must not be large enough to be owed even one representative (i.e., for j=1).")
    print("This means the condition must be false for j=1, leading to the inequality:")
    print(f"l < 1 * (s2 / k)")
    print(f"Substituting the values: {l} < s2 / {k}")
    s2 = l * k + 1
    print(f"From this, we get s2 > {l*k}. The smallest integer value for s2 is {s2}.")
    print("-" * 20)

    # Return the result as a pair (s1, s2).
    result = (s1, s2)
    print(f"The final result is the pair (s1, s2) = {result}.")

    # Final answer in the required format
    print(f"<<<{result}>>>")

if __name__ == "__main__":
    solve_election_proportionality()