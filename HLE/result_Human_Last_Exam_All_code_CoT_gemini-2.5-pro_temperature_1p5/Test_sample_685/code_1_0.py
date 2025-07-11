def f(n, m):
    """
    This function solves the 2D-Generalized NIM problem.

    It returns 1 if the first player has a winning position with a probability
    strictly more than 50% on a random n x m binary matrix, and 0 otherwise.

    The condition for the first player to have a winning probability > 50%
    is equivalent to the number of winning positions (N-positions) being strictly
    greater than the number of losing positions (P-positions).

    Analysis of this game shows that this is true for all (n, m) pairs
    except for the case n=1 and m=1, where the number of winning and losing
    positions are equal (one each).
    """

    # The condition is n*m > 1, which is true if n or m is not 1.
    # We check the only failing case explicitly.
    if n == 1 and m == 1:
        # For n=1, m=1, there is 1 losing position and 1 winning position.
        # Probability of winning is exactly 50%, not strictly more.
        result = 0
    else:
        # For all other cases, winning positions are more numerous than losing positions.
        result = 1
    
    return result

# --- Example Usage ---
# You can change these values to test other dimensions.
n_val = 5
m_val = 5

# Calculate the result using the function f
output = f(n_val, m_val)

# Print the result.
# The following print statement displays the inputs 'n' and 'm' and the computed result.
print("n =", n_val)
print("m =", m_val)
print("f(n, m) =", output)