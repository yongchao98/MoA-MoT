def solve():
    """
    Solves the problem by determining if the specified language is decidable.

    The problem asks if a decidable language L exists such that a string w is in L 
    if and only if a Turing machine T halts on w, where T is defined to halt 
    if and only if the length of w is a perfect square.

    This means L = {w | len(w) is a perfect square}.

    A language is decidable if there is an algorithm (a decider) that always halts
    and correctly determines whether any given string w is in the language.

    We can construct such an algorithm. For any input string w:
    1. Calculate its length, n = len(w).
    2. Check if n is a perfect square.

    The following code demonstrates the decision procedure for step 2. It always
    terminates because for any n, the loop variable k will eventually become
    large enough such that k*k >= n. This proves that an always-halting algorithm
    exists, and therefore, the language L is decidable.
    """
    
    # Let's test the decision procedure for a given length n.
    # The actual string w does not matter, only its length.
    n = 36
    
    print(f"Yes, such a decidable language L exists.")
    print(f"L is the set of all strings whose length is a perfect square.")
    print(f"The following algorithm decides membership in L for any string w by checking its length, n = len(w).")
    print("-" * 30)
    print(f"Running decider for a string of length n = {n}:")

    k = 0
    while True:
        square = k * k
        # The prompt requires outputting each number in the final equation.
        # This line shows the equation being tested at each step.
        print(f"Step {k}: Checking if {k} * {k} == {n}. Result: {square}")

        if square == n:
            print(f"-> {square} is equal to {n}. The length is a perfect square.")
            print(f"The algorithm halts and ACCEPTS.")
            break
        elif square > n:
            print(f"-> {square} is greater than {n}. The length is not a perfect square.")
            print(f"The algorithm halts and REJECTS.")
            break
        else: # square < n
            print(f"-> {square} is less than {n}. Continuing search...")
            k += 1

solve()