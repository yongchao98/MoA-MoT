import math

def decider_for_L(w):
    """
    This function acts as a decider for the language L.
    L = {w | len(w) is a perfect square}.
    It returns True if w is in L, and False otherwise.
    It is guaranteed to halt for any input string w.
    """
    n = len(w)
    print(f"Input string 'w' has length n = {n}.")

    if n < 0:
        print("Length cannot be negative. String is not in L.")
        return False
    
    k = 0
    while True:
        square = k * k
        print(f"Checking k = {k}: Is {k}*{k} == {n}? Current square is {square}.")
        
        if square == n:
            print(f"Yes, {n} is a perfect square ({k}*{k}).")
            print("Conclusion: The string is in language L. The decider HALTS and accepts.")
            return True
        
        if square > n:
            print(f"No, {k}*{k} ({square}) is greater than {n}. {n} cannot be a perfect square.")
            print("Conclusion: The string is not in language L. The decider HALTS and rejects.")
            return False
        
        k += 1

def main():
    """
    Main function to demonstrate the decider with example strings.
    """
    # Example 1: A string whose length is a perfect square (9)
    w1 = "101101101" # length is 9
    print("--- Running decider for w1 ---")
    decider_for_L(w1)
    
    print("\n" + "="*40 + "\n")
    
    # Example 2: A string whose length is not a perfect square (10)
    w2 = "0110101011" # length is 10
    print("--- Running decider for w2 ---")
    decider_for_L(w2)
    
    print("\n" + "="*40 + "\n")

    # Example 3: An empty string (length 0, which is 0*0)
    w3 = "" # length is 0
    print("--- Running decider for w3 ---")
    decider_for_L(w3)


if __name__ == "__main__":
    main()
