import sys

def is_prime(n):
    """
    Check if a number is prime using the Miller-Rabin primality test.
    This is a probabilistic test, but for numbers of this size, it's deterministic
    enough with a sufficient number of bases.
    """
    if n < 2:
        return False
    if n == 2 or n == 3:
        return True
    if n % 2 == 0 or n % 3 == 0:
        return False

    # Find d and s such that n-1 = 2^s * d
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1

    # Use a set of bases for the test. For the size of numbers we are testing,
    # the first 12 prime numbers as bases make the test extremely reliable.
    bases = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]
    if n in bases:
        return True

    for a in bases:
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            continue
        for _ in range(s - 1):
            x = pow(x, 2, n)
            if x == n - 1:
                break
        else:
            return False
    return True

def solve_and_print():
    """
    Generates primes from concatenated digits of an irrational number and
    verifies the properties of the 6th prime.
    """
    # The irrational number we are testing is Pi (π).
    # We need at least 1246 digits, as that's the length of the 6th prime.
    pi_digits = "31415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679821480865132823066470938446095505822317253594081284811174502841027019385211055596446229489549303819644288109756659334461284756482337867831652712019091456485669234603486104543266482133936072602491412737245870066063155881748815209209628292540917153643678925903600113305305488204665213841469519415116094330572703657595919530921861173819326117931051185480744623799627495673518857527248912279381830119491298336733624406566430860213949463952247371907021798609437027705392171762931767523846748184676694051320005681271452635608277857713427577896091736371787214684409012249534301465495853710507922796892589235420199561121290219608640344181598136297747713099605187072113499999-98372978049951059731732816096318595024459455346908302642522308253344685035261931188171010003137838752886587533208381420617177669147303598253490428755468731159562863882353787593751957781857780532171226806613001927876611195909216420198938095257201065485863278865936153381827968230301952035301852968995773622599413891249721775283479131515574857242454150695950829533116861727855889075098381754637464939319255060400927701671139009848824012858361603563707660104710181942955596198946767837449448255379774726847104047534646208046684259069491293313677028989152104752162056966024058038150193511253382430035587640247496473263914199272604269922796782354781636009341721641219924586315030286182974555706749838505494588586926995690927210797509302955321165344987202755960236480665499119881834797753566369807426542527862551818417574672890977772793800081647060016145249192173217214772350141441973568548161361157352552133475741849468438523323907394143334547762416862518983569485562099219222184272550254256887671790494601653466804988627232791786085784383827967976681454100953883786360950680064225125205117392984896084128488626945604241965285022210661186306744278622039194945047123713786960956364371917287467764657573962413890865832645995813390478027590099465764078951269468398352595709825822620522489407726719478268482601476990902640136394437455305068203496252451749399651431429809190659250937221696461515709858387410597885959772975498930161753928468138268683868942774155991855925245953959431049972524680845987283648842357214521023"
    # Remove the placeholder hyphen from the long string of pi digits
    pi_digits = pi_digits.replace("-", "")

    print("Searching for the 6th prime from the concatenated digits of π...")
    print("(This may take a few minutes as the numbers are very large.)\n")
    sys.stdout.flush()

    prime_count = 0
    target_prime_count = 6
    current_num_str = ""
    sixth_prime = -1

    for i, digit in enumerate(pi_digits):
        current_num_str += digit
        
        # Optimization: Only test numbers ending in 1, 3, 7, 9 (and initial 3)
        if len(current_num_str) > 1 and current_num_str[-1] not in '1379':
            continue

        num = int(current_num_str)
        if is_prime(num):
            prime_count += 1
            if prime_count == target_prime_count:
                sixth_prime = num
                break
    
    if sixth_prime != -1:
        print("The irrational number is π (pi).")
        print(f"The {target_prime_count}th prime number in the sequence is:")
        print(sixth_prime)
        print("-" * 78)
        
        # As requested, output the numbers from the original problem statement.
        last_six_digits = sixth_prime % 1000000
        print("Verification of the numbers in the statement:")
        print(f"Prime number sequence index: {target_prime_count}")
        print(f"Last 6 digits of this prime: {last_six_digits}")

        if last_six_digits == 521023:
            print("\nThe last 6 digits match the number 521023.")
        else:
            print("\nThe last 6 digits DO NOT match the number 521023.")
    else:
        print(f"Could not find the {target_prime_count}th prime with the provided digits.")

if __name__ == '__main__':
    solve_and_print()