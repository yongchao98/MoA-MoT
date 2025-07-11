def solve_egg_drop_time():
    """
    Calculates the egg drop time using approximations with small integers
    and prints the calculation steps.
    """

    # --- Print explanation of the plan ---
    print("Yes, your son can calculate the time. Here is a way to do it with small numbers:")
    print("-" * 30)
    print("1. First, we need the height 'h'. The formula is h = 240 * tan(22.5).")
    print("   We can use the simple fraction 5/12 to approximate tan(22.5).")
    print("\n2. Second, we need the time 't'. The formula is t^2 = (2 * h) / g.")
    print("   We can use the simple whole number 10 to approximate gravity 'g'.")
    print("\n3. Now, let's put it all together and simplify:")
    print("   t^2 = (2 * (240 * 5/12)) / 10")
    print("\nTo keep numbers small, we can think of 240 as 24 * 10:")
    print("   t^2 = (2 * (24 * 10) * 5 / 12) / 10")
    print("\nThe two 10s cancel out!")
    print("   t^2 = 2 * 24 * 5 / 12")
    print("\nSince 24 / 12 = 2, the calculation becomes much simpler:")
    print("   t^2 = 2 * 2 * 5 = 20")
    print("-" * 30)
    print("So, the time 't' is the square root of 20.")
    print("t = sqrt(20) = sqrt(4 * 5) = 2 * sqrt(5)")
    print("\n4. For the final step, we need an approximation for sqrt(5).")
    print("   A good and simple fraction for sqrt(5) is 9/4.")
    print("\nSo the final calculation for the time is:")

    # --- Final calculation ---
    num1 = 2
    num2 = 9
    den1 = 4
    result = num1 * num2 / den1
    
    # "Remember in the final code you still need to output each number in the final equation!"
    print(f"Time = {num1} * {num2} / {den1} = {result} seconds")
    print("-" * 30)

solve_egg_drop_time()