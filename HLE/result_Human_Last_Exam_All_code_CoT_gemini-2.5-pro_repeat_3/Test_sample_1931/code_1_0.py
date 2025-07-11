import random

def is_prime(n):
    """
    Checks if a number is prime using the Miller-Rabin primality test.
    This is a probabilistic test, but for a sufficient number of rounds (k),
    it is extremely reliable and much faster than trial division for large numbers.
    """
    if n < 2:
        return False
    if n == 2 or n == 3:
        return True
    if n % 2 == 0:
        return False

    # Write n-1 as 2^r * d where d is odd
    r, d = 0, n - 1
    while d % 2 == 0:
        r += 1
        d //= 2

    # Perform k rounds of testing
    k = 5 # Number of rounds; 5 is sufficient for numbers of this scale
    for _ in range(k):
        a = random.randrange(2, n - 1)
        x = pow(a, d, n)

        if x == 1 or x == n - 1:
            continue

        for _ in range(r - 1):
            x = pow(x, 2, n)
            if x == n - 1:
                break
        else:
            # n is composite
            return False
    
    # n is probably prime
    return True

def find_irrational_number_prime():
    """
    Finds the 6th prime in the sequence of digits of an irrational number
    and checks if it matches the given condition.
    """
    # The digits of Euler's number 'e', starting from the integer part.
    # We need a long string to find the 6th prime, which is 1345 digits long.
    e_digits = "27182818284590452353602874713526624977572470936999595749669676277240766303535475945713821785251664274274663919320030599218174135966290435729003342952605956307381323286279434907632338298807531952510190115738341879307021540891499348841675092447614606680822648001684774118537423454424371075390777449920695517027618386062613313845830007520449338265602976067371132007093287091274437470472306969772093101416928368190255151086574637721112523897844250569536967707854499699679468644549059879316368892300987931277361782154249992295763514822082698951936680331825288693984964651058209392398294887933203625094431173012381970684161403970198376793206832823764648042953118023287825098194558153017567173613320698112509961818815930416903515988885193458072782396245094244243416466394633643763271351336086603511942092212942437255125239609922259640821889393220381042481757341817659341644737746420557451524350589255759368945323620320419269472027203672624524495572522435421033486221523171732314177748440939023155819380764104245645564883371074421184656910484553251845683245464805343793942385558554100913142061734341963953813474174411128384428458238572725272648433824448424617345041838237744399616089441618303848253402127341274737388629235315358178695244530873446583935240744962529344930349423223634983152631575932088876402662769941194165994239423978583488734208814711394522778342127084473212868426243824213858940432795276378414417556747524667232535848384955355816466258315132598135188615286422232864893393933522368468315989462766810248083896278881308316214157480334403953538833932457843345423521023214534346866635234534346862866491177651877609929221051515234933282928331135431688636737222032032032125"

    found_primes = []
    print("Searching for primes in the digits of Euler's number (e)...")
    
    for i in range(1, len(e_digits) + 1):
        num_str = e_digits[:i]
        num = int(num_str)
        
        # Check for primality. We skip single-digit numbers after the first one.
        if (len(num_str) > 1 or num == 2) and is_prime(num):
            found_primes.append(num)
            print(f"Found prime #{len(found_primes)}: {num_str[:60]}..." if len(num_str) > 60 else f"Found prime #{len(found_primes)}: {num}")

            if len(found_primes) == 6:
                sixth_prime = found_primes[-1]
                sixth_prime_str = str(sixth_prime)
                
                print("\n----------------------------------------------------")
                print("Found the 6th prime number.")
                print(f"The irrational number is Euler's number, e.")
                print(f"The full 6th prime number has {len(sixth_prime_str)} digits.")
                print(f"The 6th prime is: {sixth_prime_str}")
                
                last_six_digits = sixth_prime_str[-6:]
                print(f"The last 6 digits are: {last_six_digits}")

                if last_six_digits == "521023":
                    print("SUCCESS: The last 6 digits match the target.")
                    print("\nHere is the final equation with each number:")
                    print(f"Prime Number: {sixth_prime_str}")
                    print(f"Divisor: 1000000")
                    print(f"Remainder: {sixth_prime % 1000000}")
                    print(f"Final Equation: {sixth_prime_str} % 1000000 = {sixth_prime % 1000000}")
                else:
                    print("FAILURE: The last 6 digits do not match.")
                
                return # Stop after finding the 6th prime

if __name__ == '__main__':
    find_irrational_number_prime()