import fractions

def solve_cathedral_riddle():
    """
    Solves the word problem from the poem "The Cathedral's Echo".
    """
    # Step 1: Identify the constant value from the poem.
    in_tune_pipes = 200
    print(f"Number of pipes that still sing pure: {in_tune_pipes}")
    print("-" * 30)

    # Step 2: Determine the fractions that make up the out-of-tune pipes.
    # The most consistent interpretation is that two groups make up the out-of-tune pipes:
    # - "One-third of pipes fell out of tune"
    # - "One-fourth descended minor scales"
    # The "two-fifths" is treated as a red herring, as it leads to contradictions with
    # the constraint that the number of out-of-tune pipes must be divisible by 28 (lcm of 7 and 4).
    fraction1 = fractions.Fraction(1, 3)
    fraction2 = fractions.Fraction(1, 4)
    total_fraction_out_of_tune = fraction1 + fraction2
    
    print(f"Fraction of pipes in the first out-of-tune group: {fraction1.numerator}/{fraction1.denominator}")
    print(f"Fraction of pipes in the second out-of-tune group: {fraction2.numerator}/{fraction2.denominator}")
    print(f"Total fraction of out-of-tune pipes (k): {total_fraction_out_of_tune.numerator}/{total_fraction_out_of_tune.denominator}")
    print("-" * 30)

    # Step 3 & 4: Formulate and solve the equation for 'O' (total out-of-tune pipes).
    # The total number of pipes, T, can be expressed as T = O + in_tune_pipes.
    # We also have the relation O = k * T.
    # Substituting T, we get the equation: O = k * (O + in_tune_pipes)
    # O = (7/12) * (O + 200)
    # 12*O = 7*(O + 200)
    # 12*O = 7*O + 1400
    # 5*O = 1400
    # O = 1400 / 5
    
    # We can solve for O directly: O * (1 - k) = k * I => O = (k * I) / (1-k)
    k = total_fraction_out_of_tune
    out_of_tune_pipes_numerator = k.numerator * in_tune_pipes
    out_of_tune_pipes_denominator = (1 - k).numerator 
    # The equation being solved is (1-k)*O = k*I, so 5/12 * O = 7/12 * 200 => 5 * O = 7 * 200
    
    equation_lhs_O_multiplier = (k.denominator - k.numerator)
    equation_rhs = k.numerator * in_tune_pipes

    print("To find the number of out-of-tune pipes (O), we solve the equation:")
    print(f"O = ({k.numerator}/{k.denominator}) * (O + {in_tune_pipes})")
    print("Which simplifies to:")
    print(f"{equation_lhs_O_multiplier} * O = {equation_rhs}")

    out_of_tune_pipes = equation_rhs / equation_lhs_O_multiplier
    print(f"\nTotal number of out-of-tune pipes (O): {int(out_of_tune_pipes)}")
    print("-" * 30)
    
    # Step 5: Calculate the final answer.
    # The tuner needs to find half of the lost pipes.
    pipes_to_find = out_of_tune_pipes / 2
    
    print("The question asks how many pipes the tuner must find when half realign.")
    print(f"This is half of the out-of-tune pipes: {int(out_of_tune_pipes)} / 2")
    print(f"\nFinal Answer: The tuner must find {int(pipes_to_find)} pipes.")

solve_cathedral_riddle()
<<<140>>>