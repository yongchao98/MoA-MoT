import requests

# The URL found in the search results for OEIS A020449
url = "https://oeis.org/A020449/b020449.txt"
response = requests.get(url)
data = response.text

lines = data.strip().split('\n')
primes_with_111 = 0
primes_without_111 = 0
total_primes = 0

for line in lines:
    # Skip comments
    if line.startswith('#'):
        continue
    
    parts = line.split()
    if len(parts) == 2:
        total_primes += 1
        prime_str = parts[1]
        if "111" in prime_str:
            primes_with_111 += 1
        else:
            primes_without_111 += 1

# The file contains exactly 1000 primes.
# So total_primes should be 1000.
# The number that will die out are the ones without "111".

# print(f"Total primes found: {total_primes}")
# print(f"Primes with '111': {primes_with_111}")
# print(f"Primes without '111': {primes_without_111}")