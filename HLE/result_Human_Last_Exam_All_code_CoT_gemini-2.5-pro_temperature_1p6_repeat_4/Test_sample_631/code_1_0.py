def solve_prime_letter_count():
    """
    This function calculates the number of 'i's in the word representations
    of the first thirteen prime numbers.
    """
    # The first 13 prime numbers and their English word representations
    primes_as_words = {
        2: "two",
        3: "three",
        5: "five",
        7: "seven",
        11: "eleven",
        13: "thirteen",
        17: "seventeen",
        19: "nineteen",
        23: "twenty-three",
        29: "twenty-nine",
        31: "thirty-one",
        37: "thirty-seven",
        41: "forty-one"
    }

    total_i_count = 0
    print("Calculating the number of 'i's in the first 13 prime numbers written as words:\n")

    # The equation part: iterate and print the count for each number
    for number, word in primes_as_words.items():
        count = word.count('i')
        print(f"{number} ('{word}') contains {count} 'i's")
        total_i_count += count

    print("\n------------------------------------")
    print(f"The final sum is the answer to the question.")
    print(f"Total number of 'i's: {total_i_count}")
    print("------------------------------------")

solve_prime_letter_count()
<<<6>>>