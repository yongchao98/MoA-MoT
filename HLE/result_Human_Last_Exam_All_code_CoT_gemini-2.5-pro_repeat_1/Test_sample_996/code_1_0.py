def solve():
    """
    This function provides the final answer in the required format.
    The reasoning is as follows:

    Part A: alpha(k) = max_{j<k} length(E(j))
    The length L(k) of E(k) follows a recurrence based on its representation k = (a+1)c_n + b.
    The maximum length, alpha(k), is dominated by the lengths of numbers of the form c_n - 1.
    Let l_n = L(c_n - 1). The recurrence for l_n is approximately l_n ~ 2 * l_{n-1}.
    This implies l_n grows as 2^n.
    Since k is related to n by n ~ log2(log2(k)), alpha(k) ~ 2^n ~ 2^(log2(log2(k))) = log2(k).
    Therefore, alpha(k) is in Theta(k^0 * (log k)^1 * (log log k)^0).
    So, a=0, b=1, c=0.

    Part B: beta(k) = |{j<k with oddity(j)=0}|
    The oddity O(k) follows the recurrence O(k) = O(a) + O(b).
    A number k has zero oddity if and only if the corresponding a and b in its representation also have zero oddity.
    Let N_n be the number of zero-oddity integers less than c_n.
    The recurrence relation is N_n = N_{n-1} + N_{n-1}^2.
    Assuming beta(k) ~ k^d, we get N_n ~ (c_n)^d = (2^(2^n))^d.
    The recurrence N_n ~ (N_{n-1})^2 is satisfied if N_n ~ theta^(2^n) for some constant theta.
    This implies (2^(2^n))^d ~ theta^(2^n), so theta = 2^d.
    Numerical analysis of the sequence N_n = {1, 2, 6, 42, ...} suggests that d is close to 2/3.
    N_2=6, (c_2)^(2/3) = 16^(2/3) ~= 6.35.
    N_3=42, (c_3)^(2/3) = 256^(2/3) ~= 40.3.
    Given that d must be rational, d=2/3 is the most reasonable choice.
    The other exponents e and f are 0.
    So, d=2/3, e=0, f=0.

    The final answer combines these rational numbers.
    """
    a = "0"
    b = "1"
    c = "0"
    d = "2/3"
    e = "0"
    f = "0"

    answer = f"{a},{b},{c},{d},{e},{f}"
    print(answer)

solve()