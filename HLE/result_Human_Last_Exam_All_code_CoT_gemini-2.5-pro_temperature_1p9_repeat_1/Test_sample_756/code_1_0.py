import numpy as np

def solve():
    """
    This function implements the analytical method to find the maximum value of |b| + |c|.
    """
    print("Let P(x) = ax^2 + bx + c.")
    print("We want to maximize |b| + |c| subject to |P(x)| <= 1 for x in [-1, 1].")
    print("We will construct a candidate polynomial Q(x) and then scale it to satisfy the condition.")

    # We select ideal values for P(1), P(-1), P(0) that would maximize b+c.
    # This corresponds to the case b>0, c>0.
    p1 = 1
    p_minus1 = -1
    p0 = 1

    # From these values, we determine the coefficients of our candidate polynomial Q(x)
    c = p0
    b = (p1 - p_minus1) / 2
    a = p1 - b - c
    
    print(f"\nConstructing a candidate polynomial Q(x) with a={a}, b={b}, c={c}.")
    print(f"Q(x) = {a}x^2 + {b}x + {c}")

    # Now, check if |Q(x)| <= 1 for all x in [-1, 1].
    # A quadratic's extrema on an interval are at the endpoints or the vertex.
    # Vertex x-coordinate:
    xv = -b / (2 * a)

    # Evaluate Q(x) at the critical points
    q_at_1 = a * (1**2) + b * 1 + c
    q_at_minus1 = a * ((-1)**2) + b * (-1) + c
    
    critical_values = [q_at_1, q_at_minus1]

    if -1 <= xv <= 1:
        q_at_xv = a * (xv**2) + b * xv + c
        critical_values.append(q_at_xv)
        print(f"The vertex is at x = {xv}, Q({xv}) = {q_at_xv}")
    
    # Find the maximum absolute value M on the interval.
    M = max(abs(v) for v in critical_values)
    print(f"The maximum absolute value of Q(x) on [-1, 1] is M = {M}.")
    
    if M <= 1:
        print("The candidate polynomial already satisfies the condition.")
        b_final, c_final = b, c
    else:
        print(f"M > 1, so we scale the polynomial by 1/M to get a valid polynomial P(x).")
        # Scale the coefficients to get the final valid polynomial P(x)
        a_final = a / M
        b_final = b / M
        c_final = c / M
        print(f"The new coefficients are: a'={a_final:.2f}, b'={b_final:.2f}, c'={c_final:.2f}")

    # Calculate the final result
    result = abs(b_final) + abs(c_final)

    print("\nThe maximum value of |b| + |c| is the sum of the absolute values of the new coefficients:")
    print(f"|{b_final:.4f}| + |{c_final:.4f}| = {result:.4f}")

solve()