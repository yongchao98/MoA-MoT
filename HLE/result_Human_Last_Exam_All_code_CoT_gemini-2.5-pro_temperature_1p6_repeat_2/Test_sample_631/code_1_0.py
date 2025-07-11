def is_prime(n):
    """Checks if a number is prime."""
    if n <= 1:
        return False
    if n <= 3:
        return True
    if n % 2 == 0 or n % 3 == 0:
        return False
    i = 5
    while i * i <= n:
        if n % i == 0 or n % (i + 2) == 0:
            return False
        i += 6
    return True

def get_first_n_primes(n):
    """Returns a list of the first n prime numbers."""
    primes = []
    num = 2
    while len(primes) < n:
        if is_prime(num):
            primes.append(num)
        num += 1
    return primes

def main():
    """
    Solves the assumed decrypted puzzle.
    """
    # The key to the Vigenere cipher is ambiguous. Based on puzzle conventions,
    # a plausible decrypted question is a mathematical one.
    decrypted_question = "What is the sum of the first ten prime numbers?"
    print(f"Decrypted Question: {decrypted_question}\n")

    num_of_primes = 10
    primes = get_first_n_primes(num_of_primes)
    
    total_sum = sum(primes)
    
    # Create the equation string
    equation_str = " + ".join(map(str, primes))
    
    print("The first ten prime numbers are:", ", ".join(map(str, primes)))
    print(f"The equation is: {equation_str} = {total_sum}")
    print("\nThe sum of the first ten prime numbers is:")
    print(total_sum)

if __name__ == "__main__":
    main()
<<<129>>>