def solve():
    """
    This problem asks for the dimension of a vector space of functions.
    Let's analyze the structure of these functions.

    A digitary function f is defined by f(sum(A_n/10^n)) = sum(T(A)_n), where T is a shortsighted map.
    A shortsighted map T has the property that T(A)_n depends only on A_n, A_{n+1}, A_{n+2}, and n.
    Let's denote this dependency as a function g_n: T(A)_n = g_n(A_n, A_{n+1}, A_{n+2}).

    The set of all sequences of functions (g_n) forms a vector space of uncountably infinite dimension, let's call it R.
    However, for f to be a well-defined function on [0, 10], its value must be the same for different decimal representations of the same number. For example, 1.0 = 0.999...
    This imposes a countable set of linear constraints on the functions g_n.

    The constraints have a specific structure: they define g_{m-2} in terms of g_k for k >= m-1. This is a backward recurrence relation.
    This means we can freely choose the functions g_n for n greater than some N, and the functions g_{N-1}, g_{N-2}, ..., g_0 are determined by these choices.

    The space of possible choices for the "tail" sequence (g_n for n >= N) is itself a vector space of uncountably infinite dimension (R). For example, we can choose the tail from a space of sequences that decay fast enough for the sum to converge. Even with such restrictions, the dimension is still uncountably infinite.

    Since we can choose the tail of the sequence (g_n) from a space of uncountable dimension, and each such choice (modulo a small-dimensional kernel) generates a valid digitary function, the resulting vector space of digitary functions must also have uncountable dimension.

    The dimension is therefore R (uncountably infinite).
    """
    # The reasoning above leads to the conclusion that the dimension is uncountably infinite.
    # The problem asks for N, R, or a natural number.
    # R represents uncountably infinite.
    answer = 'R'
    print(answer)

solve()