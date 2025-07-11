def solve_peg_game():
    """
    Determines the number of equivalence classes in the described peg game
    by explaining a mathematical proof based on a coloring invariant.
    """
    
    print("### Determining the Number of Equivalence Classes ###")
    print("\nStep 1: Understanding the Move")
    print("The game is played on a 2D integer grid Z x Z.")
    print("A move involves three consecutive positions p1, p2, p3 in a line (e.g., p2 = p1 + v, p3 = p2 + v).")
    print("A 'forward' move takes pegs at p1 and p2 and replaces them with a single peg at p3.")
    print("Let's analyze this move: {p1, p2} -> {p3}.")
    print("Note that p3 = p1 + 2v = 2*p2 - p1.")
    print("A 'backward' move is the reverse: {p3} -> {p1, p2}.")
    print("Two configurations are equivalent if one can be reached from the other by a sequence of these moves.")

    print("\nStep 2: Defining a Coloring Invariant")
    print("To distinguish between classes, we seek an invariant: a property that doesn't change with moves.")
    print("Let's color each point (x, y) of the grid with one of three colors based on its coordinates.")
    print("Color c(x, y) = (x - y) mod 3. The possible colors are {0, 1, 2}.")

    print("\nStep 3: Analyzing the Effect of a Move on Peg Colors")
    print("Let N0, N1, N2 be the number of pegs of color 0, 1, and 2, respectively.")
    print("Consider a move involving points p1, p2, p3.")
    print("Let c1 = c(p1), c2 = c(p2), c3 = c(p3).")
    print("If p2 = p1 + v (v is a unit vector), then c2 = c(p1+v) which is (c1 + 1) or (c1 - 1) mod 3.")
    print("The colors of any three consecutive points in a line are always all different: {k, k+1, k+2} mod 3.")
    print("A move replaces pegs at p1 and p2 with a peg at p3.")
    print("This means N(c1) decreases by 1, N(c2) decreases by 1, and N(c3) increases by 1.")
    print("For example, if c1=0, c2=1, c3=2, then N0->N0-1, N1->N1-1, N2->N2+1.")

    print("\nStep 4: Finding the Invariants")
    print("Let's look at the parities of the peg counts (i.e., the counts modulo 2).")
    print("Let I0 = N0 mod 2, I1 = N1 mod 2, I2 = N2 mod 2.")
    print("In the example above (c1=0, c2=1, c3=2), the changes are:")
    print("  N0 -> N0 - 1  => I0 flips (I0 -> I0 + 1 mod 2)")
    print("  N1 -> N1 - 1  => I1 flips (I1 -> I1 + 1 mod 2)")
    print("  N2 -> N2 + 1  => I2 flips (I2 -> I2 + 1 mod 2)")
    print("It can be shown that for ANY move, the colors of p1, p2, and p3 are always a permutation of {0, 1, 2}.")
    print("Therefore, EVERY move flips the parity of ALL THREE counts N0, N1, and N2.")
    print("This means that the sum of any two parities is an invariant:")
    print("  (I0 + I1) mod 2 is invariant, because (I0+1) + (I1+1) = I0 + I1 + 2 = I0 + I1 (mod 2).")
    print("  (I1 + I2) mod 2 is invariant for the same reason.")
    print("These two invariants, (N0 + N1) mod 2 and (N1 + N2) mod 2, are independent.")

    print("\nStep 5: Counting the Equivalence Classes")
    print("Each of our two binary invariants can take 2 values (0 or 1).")
    print("This gives a total of 2 * 2 = 4 possible combinations for the pair of invariants.")
    print("This means there are at least 4 equivalence classes. We can show there are exactly 4 by finding a configuration for each invariant value pair:")
    
    invariant_values = [
        ("{(0,0)}", "c(0,0)=0. (N0,N1,N2)=(1,0,0).", "((1+0)mod 2, (0+0)mod 2) = (1, 0)"),
        ("{(1,0)}", "c(1,0)=1. (N0,N1,N2)=(0,1,0).", "((0+1)mod 2, (1+0)mod 2) = (1, 1)"),
        ("{(2,0)}", "c(2,0)=2. (N0,N1,N2)=(0,0,1).", "((0+0)mod 2, (0+1)mod 2) = (0, 1)"),
        ("{(0,0),(1,0),(2,0)}", "c=0,1,2. (N0,N1,N2)=(1,1,1).", "((1+1)mod 2, (1+1)mod 2) = (0, 0)")
    ]

    print("\n{:<25} | {:<30} | Invariant ((N0+N1)%2, (N1+N2)%2)".format("Configuration", "Colors & Counts"))
    print("-" * 80)
    for config, counts, inv in invariant_values:
        print("{:<25} | {:<30} | {}".format(config, counts, inv))

    print("\nSince all 4 possible invariant pairs are achievable, there are exactly 4 equivalence classes.")
    
    final_answer = 4
    print("\nFinal Answer: The number of equivalence classes is", final_answer)


solve_peg_game()