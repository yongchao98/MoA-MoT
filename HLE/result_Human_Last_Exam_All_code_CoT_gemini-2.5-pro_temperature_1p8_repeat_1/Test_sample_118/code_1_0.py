import math

def combinations(n, k):
    """Calculates the number of combinations 'n choose k'."""
    if k < 0 or k > n:
        return 0
    # Use math.comb for efficiency and clarity (available in Python 3.8+)
    try:
        return math.comb(n, k)
    except AttributeError:
        # Fallback for older Python versions
        if k == 0 or k == n:
            return 1
        if k > n // 2:
            k = n - k
        
        res = 1
        for i in range(k):
            res = res * (n - i) // (i + 1)
        return res

def solve_coin_problem():
    """
    Calculates the maximum number of real coins we can guarantee to identify.
    """
    total_coins = 1000
    fake_coins = 4
    num_weighings = 2
    num_outcomes = 3**num_weighings

    print("Step 1: Calculate the total number of possible ways the fake coins can be distributed.")
    total_possibilities = combinations(total_coins, fake_coins)
    print(f"C({total_coins}, {fake_coins}) = {total_possibilities}")
    print("-" * 20)

    print("Step 2: Calculate the minimum number of possibilities that must be handled in the worst-case outcome.")
    print("This is the total possibilities divided by the number of outcomes from the weighings.")
    worst_case_load = math.ceil(total_possibilities / num_outcomes)
    print(f"{total_possibilities} / {num_outcomes} = {total_possibilities / num_outcomes:.2f}")
    print(f"So, the worst-case outcome must account for at least {worst_case_load} possibilities.")
    print("-" * 20)
    
    print("Step 3: Find the smallest suspect group size 'm' that can contain at least that many possibilities.")
    print("We need to find the smallest 'm' such that C(m, 4) >= " + str(worst_case_load))
    
    # We can search for m. Start with a reasonable guess based on m^4 ~ C(m,4)*24
    # m_guess = (worst_case_load * 24)**0.25 -> (4.6e9 * 24)**0.25 ~ (1.1e11)**0.25 ~ 575
    m = 570
    while True:
        current_combinations = combinations(m, fake_coins)
        if current_combinations >= worst_case_load:
            break
        m += 1
    
    min_suspects = m
    c_m_4 = combinations(min_suspects, fake_coins)
    print(f"The smallest integer 'm' is {min_suspects}.")
    print(f"Let's check: C({min_suspects - 1}, {fake_coins}) = {combinations(min_suspects - 1, fake_coins)}")
    print(f"And C({min_suspects}, {fake_coins}) = {c_m_4}, which is >= {worst_case_load}")
    print("-" * 20)

    print("Step 4: Calculate the guaranteed number of real coins.")
    print("This is the total number of coins minus the size of the worst-case suspect group 'm'.")
    guaranteed_real_coins = total_coins - min_suspects
    print(f"Guaranteed real coins = {total_coins} - {min_suspects} = {guaranteed_real_coins}")
    
    return guaranteed_real_coins

if __name__ == "__main__":
    final_answer = solve_coin_problem()
    # The final answer is wrapped separately as requested by the user prompt.
    # print(f"\nFinal Answer: {final_answer}")
    
solve_coin_problem()