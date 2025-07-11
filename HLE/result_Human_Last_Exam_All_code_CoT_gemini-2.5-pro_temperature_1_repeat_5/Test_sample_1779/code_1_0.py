import math

def check_halting_condition(input_string: str):
    """
    This function simulates a decider for the language L where strings
    must have a length that is a perfect square.

    Args:
        input_string: The string w to be checked.
    """
    length = len(input_string)
    print(f"Analyzing input string w = \"{input_string}\"...")
    print(f"Length of w is |w| = {length}")
    print("-" * 20)

    # The language L is decidable because we can create an algorithm that always
    # halts and gives a correct yes/no answer.
    # The following algorithm decides L.

    # Step 1: Handle the trivial case of length 0.
    if length == 0:
        print("Result: T halts on w.")
        print("The final equation is: 0 * 0 = 0")
        return

    # Step 2: Iterate to find an integer square root.
    # We only need to check integers i from 1 up to length.
    # The loop is guaranteed to terminate because i will eventually exceed length.
    # We can optimize by stopping when i*i > length.
    found = False
    for i in range(1, length + 1):
        square = i * i
        # If square equals length, we found the integer root. T halts.
        if square == length:
            print("Result: T halts on w.")
            print("The final equation is:")
            # Output each number in the final equation
            print(f"{i} * {i} = {length}")
            found = True
            break
        # If square exceeds length, we know length is not a perfect square. T does not halt.
        elif square > length:
            print(f"Checked up to i={i}, where {i}*{i}={square}, which is greater than {length}.")
            break

    if not found:
        print("Result: T does not halt on w, because its length is not a perfect square.")


# --- Demonstration with two examples ---

# Example 1: A string whose length is a perfect square (9)
w1 = "101101101"
check_halting_condition(w1)

print("\n" + "="*40 + "\n")

# Example 2: A string whose length is not a perfect square (10)
w2 = "0110101001"
check_halting_condition(w2)
