def solve_coin_puzzle():
    """
    This script explains and calculates the maximum number of real coins
    that can be guaranteed to be identified in two weighings.

    Problem Setup:
    - 1000 total coins
    - 4 fake coins (lighter)
    - 996 real coins
    - 1 balance scale
    - 2 weighings

    Strategy:
    1. Divide the 1000 coins into four groups:
       - Group A: 333 coins
       - Group B: 333 coins
       - Group C: 333 coins
       - Group D: 1 coin (the remainder)

    2. Perform two weighings:
       - Weighing 1: Group A vs. Group B
       - Weighing 2: Group A vs. Group C

    3. Analyze the specific outcome path that provides a guarantee:
       - Outcome 1: Weighing 1 balances (A = B).
         This means num_fake(A) = num_fake(B). Let's call this 'x'.
       - Outcome 2: Group A is lighter than Group C (A < C).
         This means num_fake(A) > num_fake(C). Let's call num_fake(C) 'y'. So, x > y.

    4. Logical Deduction:
       The total number of fakes is 4. Let num_fake(D) be 'z'.
       The equation for the fakes is:
       num_fake(A) + num_fake(B) + num_fake(C) + num_fake(D) = 4
       x + x + y + z = 4
       2x + y + z = 4

       We solve this equation with the constraints:
       - x, y, z are non-negative integers.
       - x > y (from weighing 2).
       - z <= 1 (since Group D has only one coin).

       The only integer solution that fits all constraints is x=2, y=0, z=0.

    5. Conclusion:
       This outcome proves that num_fake(C)=0 and num_fake(D)=0.
       Therefore, we can guarantee that all coins in Group C and Group D are real.

    Calculation:
    """
    group_size_n = 1000 // 3
    group_d_size = 1000 % 3
    guaranteed_real_coins = group_size_n + group_d_size

    print("The weighing strategy is:")
    print("1. Divide coins into Group A (333), B (333), C (333), and D (1).")
    print("2. Weigh A vs. B, then A vs. C.")
    print("\nIf A and B balance, and A is lighter than C, we have a guarantee.")
    print("This outcome proves that all fake coins are in groups A and B.")
    print("Therefore, groups C and D are guaranteed to be real.")
    print("\nMaximum number of guaranteed real coins is the size of Group C plus Group D.")
    print(f"The equation is: {group_size_n} + {group_d_size} = {guaranteed_real_coins}")


solve_coin_puzzle()
<<<334>>>